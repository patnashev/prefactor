#define NET_PREFACTOR_VERSION "0.8.0"
#define _SILENCE_CXX17_ALLOCATOR_VOID_DEPRECATION_WARNING

#include <stdio.h>
#include <stdlib.h>
#include "gwnum.h"
#include "cpuid.h"

#include "net.h"
#include "fermat.h"
#include "md5.h"
#include "task.h"
#include "exception.h"
#include "stage1.h"
#include "stage2.h"

void get_edecm_stage1_params(int B1, int maxSize, int *W);
void get_stage2_params(int B1, int B2, int maxSize, int *D, int *A, int *L, double *pairing);
int get_stage1_cost(int B1, int stage1param);
int get_stage2_cost(int B1, int B2, int D, int A, int L, double pairing);
int get_edecm_stage2_cost(int B1, int B2, int D, int L, int LR, double pairing);

using namespace restc_cpp;
using namespace arithmetic;

Writer* NetLogging::LoggingNetFile::get_writer()
{
    std::lock_guard<std::mutex> lock(net().upload_mutex());
    net().upload_cancel(this);
    return new Writer(std::move(_buffer));
}

void NetLogging::LoggingNetFile::on_upload()
{
    _buffer.swap(net().buffer());
    _buffer.clear();
}

void NetLogging::report(const std::string& message, int level)
{
    Logging::report(message, level);
    if (_net_level > level)
        return;
    std::unique_ptr<Writer> writer(_file.get_writer());
    writer->write(message.data(), (int)message.size());
    _file.commit_writer(*writer);
}

void NetLogging::report_factor(InputNum& input, const Giant& f)
{
    Logging::report_factor(input, f);
    _net.task()->factor = f;
}

void NetFile::on_upload()
{
    _uploading = true;
}

void NetFile::on_uploaded()
{
    _uploading = false;
}

Writer* NetFile::get_writer()
{
    std::lock_guard<std::mutex> lock(_net_ctx.upload_mutex());
    _net_ctx.upload_cancel(this);
    if (_uploading)
    {
        _buffer.swap(_net_ctx.buffer());
        _uploading = false;
    }
    return File::get_writer();
}

Reader* NetFile::get_reader()
{
    if (_buffer.empty())
    {
        boost::optional<std::string> md5;
        std::string str_data;

        // Run our example in a lambda co-routine
        auto done = _net_ctx.client()->ProcessWithPromiseT<bool>([&](Context& ctx) {
            // This is the co-routine, running in a worker-thread

            while (true)
                try
                {
                    // Construct a request to the server
                    auto reply = RequestBuilder(ctx)
                        .Get(_net_ctx.url() + "pf/" + _net_ctx.task_id() + "/" + filename())
                        .Argument("workerID", _net_ctx.worker_id())

                        // Send the request
                        .Execute();

                    md5 = reply->GetHeader("MD5");
                    str_data = reply->GetBodyAsString();

                    return true;
                }
                catch (const HttpNotFoundException&) {
                    //clog << "No file." << endl;
                    str_data.clear();
                    return true;
                }
                catch (const HttpForbiddenException&) {
                    //clog << "No task." << endl;
                    Task::abort();
                    return false;
                }
                catch (const std::exception& ex) {
                    std::clog << "File " << filename() << " download failed: " << ex.what() << std::endl;
                    ctx.Sleep(boost::posix_time::microseconds(15000000));
                    continue;
                }

            return true;
            });

        if (!done.get())
            return nullptr;

        if (hash && !str_data.empty() && md5)
        {
            char md5hash[33];
            md5_raw_input(md5hash, (unsigned char*)str_data.data(), (int)str_data.size());
            _md5hash = md5hash;
            if (md5.get() != _md5hash)
            {
                clear();
                return nullptr;
            }
        }

        _buffer.insert(_buffer.end(), str_data.begin(), str_data.end());
    }

    if (_buffer.size() < 8)
        return nullptr;
    if (*(uint32_t*)_buffer.data() != MAGIC_NUM)
        return nullptr;
    if (_buffer[4] != appid)
        return nullptr;

    return new Reader(_buffer[5], _buffer[6], _buffer[7], _buffer.data(), (int)_buffer.size(), 8);
}

void NetFile::commit_writer(Writer& writer)
{
    std::lock_guard<std::mutex> lock(_net_ctx.upload_mutex());
    if (_uploading)
    {
        _buffer.swap(_net_ctx.buffer());
        _uploading = false;
    }
    if (hash)
        _md5hash = writer.hash_str();
    _buffer = std::move(writer.buffer());
    _net_ctx.upload(this);
}

void NetFile::clear()
{
    std::lock_guard<std::mutex> lock(_net_ctx.upload_mutex());
    _net_ctx.upload_cancel(this);
    if (_uploading)
    {
        _buffer.swap(_net_ctx.buffer());
        _uploading = false;
    }
    _buffer.clear();
    _md5hash.clear();
    _net_ctx.upload(this);
}

void NetContext::upload_cancel(NetFile* file)
{
    auto it = _upload_queue.begin();
    while (it != _upload_queue.end())
        if (*it == file)
            it = _upload_queue.erase(it);
        else
            it++;
}

void NetContext::upload_wait()
{
    if (_uploadF.valid())
        _uploadF.get();
}

void NetContext::upload(NetFile* file)
{
    _upload_queue.push_back(file);
    if (_uploadF.valid() || _task->aborted)
        return;

    _uploadF = _putter->ProcessWithPromise([this](Context& ctx) {

        NetFile* file = nullptr;
        std::string put_url;
        std::string md5;
        char* data;
        size_t size;
        while (true)
        {
            if (file == nullptr)
            {
                std::lock_guard<std::mutex> lock(_upload_mutex);
                if (_upload_queue.empty() || _task->aborted)
                {
                    _uploadF = std::future<void>();
                    return;
                }
                file = _upload_queue.front();
                _upload_queue.pop_front();
                put_url = url() + "pf/" + task_id() + "/" + file->filename();
                data = file->buffer().data();
                size = file->buffer().size();
                md5 = file->md5hash();
                file->on_upload();
            }

            try
            {
                auto properties = _putter->GetConnectionProperties();
                properties->afterWriteFn = [&]()
                {
                    std::lock_guard<std::mutex> lock(_upload_mutex);
                    file->on_uploaded();
                };

                RequestBuilder(ctx)
                    .Put(put_url)
                    .Argument("md5", md5)
                    .Argument("workerID", worker_id())
                    .Argument("uptime", uptime())
                    .Argument("fft_desc", _task->fft_desc)
                    .Argument("fft_len", _task->fft_len)
                    .Argument("progress", _task->progress)
                    .Argument("timer", _task->timer)
                    .Body(std::unique_ptr<RequestBody>(new RequestBodyData(data, size)))
                    .Execute();

                properties->afterWriteFn = nullptr;
                file = nullptr;
            }
            catch (const HttpAuthenticationException&) {
                std::clog << "Task timed out." << std::endl;
                _task->aborted = true;
                Task::abort();
                file = nullptr;
            }
            catch (const HttpForbiddenException&) {
                std::clog << "Task not found." << std::endl;
                _task->aborted = true;
                Task::abort();
                file = nullptr;
            }
            catch (const std::exception& ex) {
                std::clog << "Upload to " << put_url << " failed: " << ex.what() << std::endl;
                ctx.Sleep(boost::posix_time::microseconds(15000000));
            }
        }
    });
}

void NetContext::done()
{
    _putter->CloseWhenReady(true);
}


/*
extern "C" void net_report_init(const char *fft_desc, int fft_len, nethandle net)
{
	if (net != NULL)
	{
		NetContext *ctx = (NetContext*)net;
		ctx->task->fft_desc = fft_desc;
		ctx->task->fft_len = fft_len;
        return;
	}
}

extern "C" void net_report_progress(int progress, int timer, nethandle net)
{
    if (net != NULL)
    {
        NetContext *ctx = (NetContext*)net;
        ctx->task->progress = progress;
        ctx->task->timer = timer;
        return;
    }
}
*/
int net_main(int argc, char *argv[])
{
    int i;
    GWState gwstate;
    std::string url;
    std::string worker_id;
    int log_level = Logging::LEVEL_INFO;
    int net_log_level = Logging::LEVEL_WARNING;
    int maxMem = 2048;

    for (i = 1; i < argc; i++)
        if (argv[i][0] == '-' && argv[i][1])
        {
            switch (argv[i][1])
            {
            case 't':
                if (argv[i][2] && isdigit(argv[i][2]))
                    gwstate.thread_count = atoi(argv[i] + 2);
                else if (!argv[i][2] && i < argc - 1)
                {
                    i++;
                    gwstate.thread_count = atoi(argv[i]);
                }
                else
                    break;
                if (gwstate.thread_count == 0 || gwstate.thread_count > 64)
                    gwstate.thread_count = 1;
                continue;

            case 'M':
                if (argv[i][2] && isdigit(argv[i][2]))
                    maxMem = atoi(argv[i] + 2);
                else if (!argv[i][2] && i < argc - 1)
                {
                    i++;
                    maxMem = atoi(argv[i]);
                }
                else
                    break;
                continue;

            case 'i':
                if (argv[i][2])
                    worker_id = argv[i] + 2;
                else if (!argv[i][2] && i < argc - 1)
                {
                    i++;
                    worker_id = argv[i];
                }
                else
                    break;
                continue;
            }

            if (strcmp(argv[i], "-v") == 0)
            {
                printf("Net-Prefactor version " NET_PREFACTOR_VERSION ", Gwnum library version " GWNUM_VERSION "\n");
                return 0;
            }
        }
        else
        {
            url = argv[i];
        }
    if (url.empty() || url.find("http://") != 0 || worker_id.empty())
    {
        printf("Usage: prefactor -net -i WorkerID http://host:port/api/\n");
        return 0;
    }

    NetContext net(url, worker_id, log_level, net_log_level);
    Logging& logging = net.logging();

	// Set the log-level to a reasonable value
	boost::log::core::get()->set_filter
	(
#ifdef _DEBUG
		boost::log::trivial::severity >= boost::log::trivial::warning
#else
		boost::log::trivial::severity >= boost::log::trivial::error
#endif // DEBUG
	);

	while (true)
	{
        if (Task::abort_flag())
            return 1;
        logging.set_prefix("");
        net.task().reset(new PFTask());

		// Run our example in a lambda co-routine
		auto done = net.client()->ProcessWithPromiseT<bool>([&](Context& ctx) {
			// This is the co-routine, running in a worker-thread

			try
			{

				// Construct a request to the server
				SerializeFromJson(*net.task(), RequestBuilder(ctx)
					.Post(net.url() + "pf/new")
					.Argument("workerID", net.worker_id())
					.Argument("uptime", net.uptime())
					.Argument("version", NET_PREFACTOR_VERSION)

					// Send the request
					.Execute()
				);
			}
			catch (const std::exception& ex) {
                std::clog << "Task acquisition failed: " << ex.what() << std::endl;
				return false;
			}

            std::cout << net.task_id() << std::endl;

			return true;
		});

		if (!done.get())
		{
			std::this_thread::sleep_for(std::chrono::minutes(1));
			continue;
		}

        InputNum input(net.task()->sk, net.task()->sb, net.task()->n, net.task()->c);
        input.setup(gwstate);
        logging.info("Using %s.\n", gwstate.fft_description.data());
        net.task()->fft_desc = gwstate.fft_description;
        net.task()->fft_len = gwstate.fft_length;
        int maxSize = (int)(maxMem/(gwnum_size(gwstate.gwdata())/1048576.0));

        logging.set_prefix(net.task_id() + ", ");
        net.task()->factor = 1;
        std::string dhash = "";
        if (!net.task()->factors.empty())
        {
            Giant factors;
            factors = net.task()->factors;
            if (*gwstate.N % factors != 0)
                logging.error("Factors do not divide the number.\n");
            *gwstate.N /= factors;
        }
        if (net.task()->b2 < net.task()->b1)
            net.task()->b2 = net.task()->b1;

        NetFile file_input(net, "input", gwstate.fingerprint);
        NetFile file_output(net, "output", gwstate.fingerprint);
        NetFile file_checkpoint(net, "checkpoint", gwstate.fingerprint);
        NetFile file_checkpoint_m1(net, "checkpoint.m1", gwstate.fingerprint);
        NetFile file_checkpoint_m12(net, "checkpoint.m12", gwstate.fingerprint);
        NetFile file_checkpoint_m2(net, "checkpoint.m2", gwstate.fingerprint);
        NetFile file_checkpoint_p1(net, "checkpoint.p1", gwstate.fingerprint);
        NetFile file_checkpoint_p12(net, "checkpoint.p12", gwstate.fingerprint);
        NetFile file_checkpoint_p2(net, "checkpoint.p2", gwstate.fingerprint);
        NetFile file_checkpoint_ed1(net, "checkpoint.ed1", gwstate.fingerprint);
        NetFile file_checkpoint_ed12(net, "checkpoint.ed12", gwstate.fingerprint);
        NetFile file_checkpoint_ed2(net, "checkpoint.ed2", gwstate.fingerprint);

        try
        {
            if (net.task()->type == "Fermat")
            {
                Fermat fermat(net.task()->n, net.task_id(), gwstate, logging);

                if (!fermat.read_points(file_input))
                {
                    logging.error("Input is missing or corrupted.\n");
                    std::this_thread::sleep_for(std::chrono::minutes(1));
                    net.task()->aborted = true;
                    throw TaskAbortException();
                }

                if (net.task()->b1 > fermat.B0())
                {
                    fermat.read_state(file_checkpoint, net.task()->b1);
                    fermat.stage1(net.task()->b1, file_checkpoint, file_output);
                }
                else
                {
                    if (net.task()->b1 != 0)
                        logging.warning("Input is at a higher B1.\n");
                }

                dhash = fermat.verify(true);
            }
            else
            {
                int B1 = (int)net.task()->b1;
                int B2 = (int)net.task()->b2;
                int D, A, L;
                double pairing;
                PrimeList primes(B1 + 100);

                if (net.task()->type == "P-1")
                {
                    get_stage2_params(B1, B2, maxSize, &D, &A, &L, &pairing);
                    logging.progress().add_stage(get_stage1_cost(B1, 0));
                    logging.progress().add_stage(get_stage2_cost(B1, B2, D, A, L, pairing));
                    PP1Stage1::State* interstate = read_state<PP1Stage1::State>(&file_checkpoint_m12);
                    if (interstate == nullptr)
                    {
                        PM1Stage1 stage1(primes, B1);
                        stage1.init(&input, &gwstate, &file_checkpoint_m1, &logging);
                        stage1.run();
                        if (!stage1.success() && B2 > B1)
                        {
                            interstate = new PP1Stage1::State();
                            interstate->V() = std::move(stage1.V());
                            file_checkpoint_m12.write(*interstate);
                        }
                    }
                    logging.progress().next_stage();
                    if (interstate != nullptr)
                    {
                        PP1Stage2 stage2(logging, primes, B1, B2, D, A, L);
                        stage2.init(&input, &gwstate, &file_checkpoint_m2, &logging, interstate->V(), true);
                        stage2.run();
                    }
                    logging.progress().next_stage();
                }
                if (net.task()->type == "P+1")
                {
                    if (net.task()->seed.empty())
                        net.task()->seed = "2/7";
                    get_stage2_params(B1, B2, maxSize, &D, &A, &L, &pairing);
                    logging.progress().add_stage(get_stage1_cost(B1, 1));
                    logging.progress().add_stage(get_stage2_cost(B1, B2, D, A, L, pairing));
                    PP1Stage1::State* interstate = read_state<PP1Stage1::State>(&file_checkpoint_p12);
                    if (interstate == nullptr)
                    {
                        PP1Stage1 stage1(primes, B1, net.task()->seed);
                        stage1.init(&input, &gwstate, &file_checkpoint_p1, &logging);
                        stage1.run();
                        if (!stage1.success() && B2 > B1)
                        {
                            interstate = new PP1Stage1::State(std::move(*stage1.state()));
                            file_checkpoint_p12.write(*interstate);
                        }
                    }
                    logging.progress().next_stage();
                    if (interstate != nullptr)
                    {
                        PP1Stage2 stage2(logging, primes, B1, B2, D, A, L);
                        stage2.init(&input, &gwstate, &file_checkpoint_p2, &logging, interstate->V(), false);
                        stage2.run();
                    }
                    logging.progress().next_stage();
                }
                if (net.task()->type == "EdECM")
                {
                    Giant X, Y, Z, T, EdD;
                    {
                        GWArithmetic gw(gwstate);
                        EdwardsArithmetic ed(gw.carefully());
                        EdPoint P(ed);
                        GWNum ed_d(gw.carefully());

                        if (net.task()->seed.empty() || net.task()->seed == "curve2x8")
                            P = ed.from_small(17, 19, 17, 33, &ed_d);
                        else if (net.task()->seed == "curve12")
                            P = ed.from_small(5, 23, -1, 7, &ed_d);
                        else
                        {
                            int curveSeed;
                            if (net.task()->seed == "random")
                            {
                                double timer = getHighResTimer();
                                curveSeed = *(int *)&timer;
                            }
                            else
                                curveSeed = std::stoi(net.task()->seed);
                            try
                            {
                                P = ed.gen_curve(curveSeed, &ed_d);
                            }
                            catch (const ArithmeticException&)
                            {
                                logging.error("Invalid curve.\n");
                                std::this_thread::sleep_for(std::chrono::minutes(1));
                                net.task()->aborted = true;
                                throw TaskAbortException();
                            }
                        }
                        Giant tmp;
                        tmp = ed.jinvariant(ed_d);
                        dhash = tmp.to_string();
                        P.serialize(X, Y, Z, T);
                        EdD = ed_d;
                    }
                    int W;
                    get_edecm_stage1_params(B1, maxSize, &W);
                    get_stage2_params(B1, B2, maxSize, &D, &A, &L, &pairing);
                    logging.progress().add_stage(get_stage1_cost(B1, W));
                    logging.progress().add_stage(get_stage2_cost(B1, B2, D, A, L, pairing));
                    EdECMStage1::State* interstate = read_state<EdECMStage1::State>(&file_checkpoint_ed12);
                    if (interstate == nullptr)
                    {
                        EdECMStage1 stage1(primes, B1, W);
                        stage1.init(&input, &gwstate, &file_checkpoint_ed1, &logging, X, Y, Z, T, EdD);
                        stage1.run();
                        if (!stage1.success() && B2 > B1)
                        {
                            interstate = new EdECMStage1::State(std::move(*stage1.state()));
                            file_checkpoint_ed12.write(*interstate);
                        }
                    }
                    logging.progress().next_stage();
                    if (interstate != nullptr)
                    {
                        EdECMStage2 stage2(logging, primes, B1, B2, 210, 5, 20);
                        stage2.init(&input, &gwstate, &file_checkpoint_ed2, &logging, interstate->X(), interstate->Y(), interstate->Z(), interstate->T(), EdD);
                        stage2.run();
                    }
                    logging.progress().next_stage();
                }
            }
        }
        catch (const NoInverseException& e)
        {
            logging.report_factor(input, e.divisor);
        }
        catch (const TaskAbortException&)
        {
        }

        gwstate.done();

        net.upload_wait();
        if (net.task()->aborted)
        {
            Task::abort_reset();
            continue;
        }
        if (Task::abort_flag())
            return 1;

		// Run our example in a lambda co-routine
		auto doneRet = net.client()->ProcessWithPromise([&](Context& ctx) {
			// This is the co-routine, running in a worker-thread

			int i;
			for (i = 0; i < 10; i++)
				try
				{
					// Construct a request to the server
					RequestBuilder(ctx)
						.Post(net.url() + "pf/res/" + net.task_id())
						.Argument("workerID", net.worker_id())
                        .Argument("result", net.task()->factor.to_string())
                        .Argument("dhash", dhash)
                        .Argument("version", NET_PREFACTOR_VERSION)

						// Send the request
						.Execute();
					break;
				}
				catch (const std::exception& ex) {
					std::clog << "Can't upload result: " << ex.what() << std::endl;
					ctx.Sleep(boost::posix_time::microseconds(i*5000000));
				}
		});

		doneRet.get();
	}

    net.done();

    return 0;
}

