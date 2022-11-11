#define NET_PREFACTOR_VERSION "0.10.0"
#define _SILENCE_CXX17_ALLOCATOR_VOID_DEPRECATION_WARNING

#include <stdio.h>
#include <stdlib.h>
#include "gwnum.h"
#include "cpuid.h"

#include "net.h"
#include "factoring.h"
#include "md5.h"
#include "task.h"
#include "exception.h"
#include "stage1.h"
#include "stage2.h"
#include "stage2poly.h"
#include "params.h"

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
    if (level >= _net_level && level != LEVEL_RESULT)
        _file.write_text(message);
    if (level == LEVEL_RESULT)
        _file_factors.write_text(_prefix + message);
}

void NetLogging::report_factor(InputNum& input, const Giant& f)
{
    Logging::report_factor(input, f);
    if (_net.task()->factor == "1")
        _net.task()->factor.clear();
    else
        _net.task()->factor += ',';
    _net.task()->factor += f.to_string();
}

void NetLogging::report_param(const std::string& name, int value)
{
    if (name == "fft_len")
        _net.task()->fft_len = value;
    if (name == "W")
        _net.task()->W = value;
    if (name == "D")
        _net.task()->D = value;
    if (name == "degree")
        _net.task()->degree = value;
}

void NetLogging::report_param(const std::string& name, const std::string& value)
{
    if (name == "fft_desc")
        _net.task()->fft_desc = value;
}

void NetLogging::progress_save()
{
    _net.task()->progress = progress().progress_total();
    _net.task()->time = progress().time_total();
    _net.task()->time_op = progress().time_op()*1000;
}

void NetFile::on_upload()
{
    _uploading = true;
}

void NetFile::on_uploaded()
{
    _uploading = false;
    if (_free_buffer)
    {
        _free_buffer = false;
        std::vector<char>().swap(_buffer);
    }
}

File* NetFile::add_child(const std::string& name, uint32_t fingerprint)
{
    _children.emplace_back(new NetFile(_net_ctx, _filename + "." + name, fingerprint));
    return _children.back().get();
}

Writer* NetFile::get_writer()
{
    std::lock_guard<std::mutex> lock(_net_ctx.upload_mutex());
    _free_buffer = false;
    _net_ctx.upload_cancel(this);
    if (_uploading)
    {
        _buffer.swap(_net_ctx.buffer());
        _uploading = false;
    }
    Writer* writer = new Writer(std::move(_buffer));
    writer->buffer().clear();
    return writer;
}

void NetFile::read_buffer()
{
    if (_free_buffer)
    {
        std::lock_guard<std::mutex> lock(_net_ctx.upload_mutex());
        _free_buffer = false;
    }
    if (!_buffer.empty())
        return;

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
        return;

    if (hash && !str_data.empty() && md5)
    {
        char md5hash[33];
        md5_raw_input(md5hash, (unsigned char*)str_data.data(), (int)str_data.size());
        _md5hash = md5hash;
        if (md5.get() != _md5hash)
        {
            clear();
            return;
        }
    }

    _buffer.insert(_buffer.end(), str_data.begin(), str_data.end());
}

void NetFile::commit_writer(Writer& writer)
{
    std::lock_guard<std::mutex> lock(_net_ctx.upload_mutex());
    _free_buffer = false;
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

void NetFile::free_buffer()
{
    std::lock_guard<std::mutex> lock(_net_ctx.upload_mutex());
    if (_net_ctx.upload_queued(this))
    {
        _free_buffer = true;
        return;
    }
    if (_uploading)
    {
        _buffer.swap(_net_ctx.buffer());
        _uploading = false;
    }
    _free_buffer = false;
    std::vector<char>().swap(_buffer);
}

void NetFile::clear(bool recursive)
{
    if (recursive)
        for (auto it = _children.begin(); it != _children.end(); it++)
            (*it)->clear(true);

    std::lock_guard<std::mutex> lock(_net_ctx.upload_mutex());
    _net_ctx.upload_cancel(this);
    if (_uploading)
    {
        _buffer.swap(_net_ctx.buffer());
        _uploading = false;
    }
    _free_buffer = false;
    std::vector<char>().swap(_buffer);
    _md5hash.clear();
    _net_ctx.upload(this);
}

bool NetContext::upload_queued(NetFile* file)
{
    for (auto it = _upload_queue.begin(); it != _upload_queue.end(); it++)
        if (*it == file)
            return true;
    return false;
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
    std::future<void> localF;
    {
        std::lock_guard<std::mutex> lock(_upload_mutex);
        localF = std::move(_uploadF);
    }
    if (localF.valid())
        localF.get();
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
                RequestBuilder(ctx)
                    .Put(put_url)
                    .Argument("md5", md5)
                    .Argument("workerID", worker_id())
                    .Argument("uptime", uptime())
                    .Argument("fft_desc", _task->fft_desc)
                    .Argument("fft_len", _task->fft_len)
                    .Argument("W", _task->W)
                    .Argument("D", _task->D)
                    .Argument("degree", _task->degree)
                    .Argument("progress", std::to_string(_task->progress))
                    .Argument("time", std::to_string(_task->time))
                    .Argument("time_op", std::to_string(_task->time_op))
                    .Body(std::unique_ptr<RequestBody>(new RequestBodyData(data, size)))
                    .Execute();

                std::lock_guard<std::mutex> lock(_upload_mutex);
                file->on_uploaded();
                file = nullptr;
            }
            catch (const HttpAuthenticationException&) {
                std::clog << "Task timed out." << std::endl;
                _task->aborted = true;
                Task::abort();
                std::lock_guard<std::mutex> lock(_upload_mutex);
                file->on_uploaded();
                file = nullptr;
            }
            catch (const HttpForbiddenException&) {
                std::clog << "Task not found." << std::endl;
                _task->aborted = true;
                Task::abort();
                std::lock_guard<std::mutex> lock(_upload_mutex);
                file->on_uploaded();
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


int net_main(int argc, char *argv[])
{
    int i;
    GWState gwstate;
    std::string url;
    std::string worker_id;
    int log_level = Logging::LEVEL_INFO;
    int net_log_level = Logging::LEVEL_WARNING;
    uint64_t maxMem = 2048*1048576ULL;
    int polyThreads = 1;
    int polyMemModel = 0;
    int disk_write_time = Task::DISK_WRITE_TIME;

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
                    maxMem = InputNum::parse_numeral(argv[i] + 2);
                else if (!argv[i][2] && i < argc - 1)
                {
                    i++;
                    maxMem = InputNum::parse_numeral(argv[i]);
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

            if (i < argc - 1 && strcmp(argv[i], "-L3") == 0)
            {
                i++;
                PolyMult::L3_CACHE_MB = atoi(argv[i]);
            }
            else if (strcmp(argv[i], "-poly") == 0)
            {
                while (true)
                    if (i < argc - 2 && strcmp(argv[i + 1], "threads") == 0)
                    {
                        i += 2;
                        polyThreads = atoi(argv[i]);
                    }
                    else if (i < argc - 1 && argv[i + 1][0] == 't')
                    {
                        i++;
                        polyThreads = atoi(argv[i] + 1);
                    }
                    else if (i < argc - 2 && strcmp(argv[i + 1], "mem") == 0)
                    {
                        i += 2;
                        if (strcmp(argv[i], "lowest") == 0)
                            polyMemModel = -2;
                        if (strcmp(argv[i], "low") == 0)
                            polyMemModel = -1;
                        if (strcmp(argv[i], "normal") == 0)
                            polyMemModel = 0;
                        if (strcmp(argv[i], "high") == 0)
                            polyMemModel = 1;
                        if (strcmp(argv[i], "highest") == 0)
                            polyMemModel = 2;
                    }
                    else
                        break;
            }
            else if (strcmp(argv[i], "-time") == 0)
            {
                while (true)
                    if (i < argc - 2 && strcmp(argv[i + 1], "write") == 0)
                    {
                        i += 2;
                        disk_write_time = atoi(argv[i]);
                    }
                    else if (i < argc - 2 && strcmp(argv[i + 1], "progress") == 0)
                    {
                        i += 2;
                        Task::PROGRESS_TIME = atoi(argv[i]);
                    }
                    else
                        break;
            }
            else if (i < argc - 1 && strcmp(argv[i], "-log") == 0)
            {
                i++;
                if (strcmp(argv[i], "debug_internal") == 0)
                    log_level = Logging::LEVEL_DEBUG - 1;
                if (strcmp(argv[i], "debug") == 0)
                    log_level = Logging::LEVEL_DEBUG;
                if (strcmp(argv[i], "info") == 0)
                    log_level = Logging::LEVEL_INFO;
                if (strcmp(argv[i], "warning") == 0)
                    log_level = Logging::LEVEL_WARNING;
                if (strcmp(argv[i], "error") == 0)
                    log_level = Logging::LEVEL_ERROR;
            }
            else if (strcmp(argv[i], "-v") == 0)
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

			return true;
		});
		if (!done.get())
		{
			std::this_thread::sleep_for(std::chrono::minutes(1));
			continue;
		}
        logging.info("%s\n", net.task_id().data());

        NetFile file_number(net, "number", 0);
        InputNum input;
        if (net.task()->n > 0)
            input.init(net.task()->sk, net.task()->sb, net.task()->n, net.task()->c);
        else if (!input.read(file_number))
        {
            logging.error("Number file is missing or corrupted.\n");
            std::this_thread::sleep_for(std::chrono::minutes(1));
            continue;
        }
        if (!net.task()->factors.empty())
            gwstate.known_factors = net.task()->factors;
        if (net.task()->options.find("FFT_Increment") != net.task()->options.end())
            gwstate.next_fft_count = std::stoi(net.task()->options["FFT_Increment"]);
        input.setup(gwstate);
        logging.info("Using %s.\n", gwstate.fft_description.data());
        net.task()->fft_desc = gwstate.fft_description;
        net.task()->fft_len = gwstate.fft_length;
        net.task()->W = net.task()->D = net.task()->degree = 0;
        int maxSize = (int)(maxMem/(gwnum_size(gwstate.gwdata())));

        logging.set_prefix(net.task_id() + ", ");
        logging.progress() = Progress();
        logging.progress().time_init(net.task()->time);
        net.task()->factor = "1";
        std::string dhash = "";
        if (net.task()->b2 < net.task()->b1)
            net.task()->b2 = net.task()->b1;
        bool poly = false;
        if (net.task()->options.find("poly") != net.task()->options.end())
            poly = net.task()->options["poly"] == "true";
        int polyCheck = 0;
        if (net.task()->options.find("poly_check") != net.task()->options.end())
            polyCheck = std::stoi(net.task()->options["poly_check"]);
        if (net.task()->options.find("write_time") != net.task()->options.end())
            Task::DISK_WRITE_TIME = std::stoi(net.task()->options["write_time"]);
        else
            Task::DISK_WRITE_TIME = disk_write_time;

        std::list<NetFile> files;
        try
        {
            if (net.task()->type == "Factoring" || net.task()->type == "Factoring2")
            {
                Factoring factoring(input, gwstate, logging);
                NetFile& file_input = files.emplace_back(net, "input", gwstate.fingerprint);
                NetFile& file_output = files.emplace_back(net, "output", gwstate.fingerprint);
                NetFile& file_dhash = files.emplace_back(net, "dhash", 0);

                if (!factoring.read_points(file_input))
                {
                    if (net.task()->options.find("generate") != net.task()->options.end())
                    {
                        std::string gen = net.task()->options["generate"];
                        int seed = std::stoi(gen);
                        int count = 1;
                        if (gen.find(",") != std::string::npos)
                            count = std::stoi(gen.substr(gen.find(",") + 1));
                        dhash = factoring.generate(seed, count);
                        factoring.write_points(file_input);
                        file_dhash.write_text(dhash);
                    }
                    else
                    {
                        logging.error("Input is missing or corrupted.\n");
                        std::this_thread::sleep_for(std::chrono::minutes(1));
                        net.task()->aborted = true;
                        throw TaskAbortException();
                    }
                }

                if (net.task()->options.find("range") != net.task()->options.end())
                {
                    std::string range = net.task()->options["range"];
                    int range_offset = std::stoi(range);
                    int range_count = 1;
                    if (range.find(",") != std::string::npos)
                        range_count = std::stoi(range.substr(range.find(",") + 1));
                    factoring.split(range_offset, range_count, factoring);
                    Factoring factoring_range(input, gwstate, logging);
                    factoring.copy(factoring_range);
                    dhash = factoring_range.verify(false);
                    file_dhash.write_text(dhash);
                }

                if (net.task()->b1 > factoring.B1())
                {
                    uint64_t B1max = net.task()->b1;
                    if (net.task()->options.find("B1max") != net.task()->options.end())
                        B1max = InputNum::parse_numeral(net.task()->options["B1max"]);

                    NetFile& file_checkpoint = files.emplace_back(net, "checkpoint", File::unique_fingerprint(gwstate.fingerprint, std::to_string(factoring.B1()) + "." + std::to_string(B1max)));
                    factoring.read_state(file_checkpoint, net.task()->b1);
                    factoring.stage1(net.task()->b1, B1max, maxMem, file_checkpoint, file_output);
                    dhash = factoring.verify(true);
                }
                else
                {
                    if (net.task()->b1 != 0 && net.task()->b2 <= factoring.B1())
                        logging.warning("Input is at a higher B1.\n");
                }

                if (net.task()->b2 > factoring.B1())
                {
                    if (!dhash.empty())
                        dhash += "-";
                    NetFile& file_results = files.emplace_back(net, "results", 0);
                    NetFile& file_checkpoint = files.emplace_back(net, "checkpoint", 0);
                    logging.progress().add_stage((int)factoring.points().size());
                    dhash += factoring.stage2(net.task()->b2, maxMem, poly, polyThreads, polyMemModel, polyCheck, &file_checkpoint, file_results);
                    logging.progress().next_stage();
                }
            }
            else
            {
                PrimeList primes((int)net.task()->b2 + 100);

                if (net.task()->type == "P-1")
                {
                    std::unique_ptr<PM1Params> params_pm1;
                    params_pm1.reset(new PM1Params(net.task()->b1, net.task()->b2, maxSize, poly, polyThreads, polyMemModel));
                    logging.progress().add_stage(params_pm1->stage1_cost());
                    logging.progress().add_stage(params_pm1->stage2_cost());

                    uint32_t fingerprint = File::unique_fingerprint(gwstate.fingerprint, std::to_string(params_pm1->B1));
                    NetFile& file1 = files.emplace_back(net, "checkpoint.m1", fingerprint);
                    NetFile& file12 = files.emplace_back(net, "checkpoint.m12", fingerprint);
                    NetFile& file2 = files.emplace_back(net, "checkpoint.m2", File::unique_fingerprint(fingerprint, std::to_string(params_pm1->B2)));
                    PP1Stage1::State* interstate = read_state<PP1Stage1::State>(&file12);
                    if (interstate == nullptr)
                    {
                        PM1Stage1 stage1((int)params_pm1->B1);
                        stage1.init(&input, &gwstate, &file1, &logging);
                        stage1.run();
                        if (!stage1.success() && params_pm1->B2 > params_pm1->B1)
                        {
                            interstate = new PP1Stage1::State();
                            interstate->V() = std::move(stage1.V());
                            file12.write(*interstate);
                        }
                    }
                    logging.progress().next_stage();
                    if (interstate != nullptr)
                    {
                        std::unique_ptr<Stage2> stage2;
                        if (params_pm1->PolyPower == 0)
                            stage2.reset(new PP1Stage2(params_pm1->B1, params_pm1->B2, params_pm1->D, params_pm1->A, params_pm1->L, logging));
                        else
                            stage2.reset(new PP1Stage2Poly(params_pm1->B1, params_pm1->B2, params_pm1->D, params_pm1->PolyPower, polyThreads, polyMemModel, polyCheck > 0));
                        dynamic_cast<IPP1Stage2*>(stage2.get())->init(&input, &gwstate, &file2, &logging, interstate->V(), true);
                        stage2->run();
                    }
                    logging.progress().next_stage();
                }
                if (net.task()->type == "P+1")
                {
                    if (net.task()->seed.empty())
                        net.task()->seed = "2/7";
                    std::unique_ptr<PP1Params> params_pp1;
                    params_pp1.reset(new PP1Params(net.task()->b1, net.task()->b2, maxSize, poly, polyThreads, polyMemModel));
                    logging.progress().add_stage(params_pp1->stage1_cost());
                    logging.progress().add_stage(params_pp1->stage2_cost());

                    uint32_t fingerprint = File::unique_fingerprint(gwstate.fingerprint, net.task()->seed + "." + std::to_string(params_pp1->B1));
                    NetFile& file1 = files.emplace_back(net, "checkpoint.p1", fingerprint);
                    NetFile& file12 = files.emplace_back(net, "checkpoint.p12", fingerprint);
                    NetFile& file2 = files.emplace_back(net, "checkpoint.p2", File::unique_fingerprint(fingerprint, std::to_string(params_pp1->B2)));
                    PP1Stage1::State* interstate = read_state<PP1Stage1::State>(&file12);
                    if (interstate == nullptr)
                    {
                        PP1Stage1 stage1((int)params_pp1->B1, net.task()->seed);
                        stage1.init(&input, &gwstate, &file1, &logging);
                        stage1.run();
                        if (!stage1.success() && params_pp1->B2 > params_pp1->B1)
                        {
                            interstate = new PP1Stage1::State(std::move(*stage1.state()));
                            file12.write(*interstate);
                        }
                    }
                    logging.progress().next_stage();
                    if (interstate != nullptr)
                    {
                        std::unique_ptr<Stage2> stage2;
                        if (params_pp1->PolyPower == 0)
                            stage2.reset(new PP1Stage2(params_pp1->B1, params_pp1->B2, params_pp1->D, params_pp1->A, params_pp1->L, logging));
                        else
                            stage2.reset(new PP1Stage2Poly(params_pp1->B1, params_pp1->B2, params_pp1->D, params_pp1->PolyPower, polyThreads, polyMemModel, polyCheck > 0));
                        dynamic_cast<IPP1Stage2*>(stage2.get())->init(&input, &gwstate, &file2, &logging, interstate->V(), false);
                        stage2->run();
                    }
                    logging.progress().next_stage();
                }
                if (net.task()->type == "EdECM")
                {
                    std::string jinvariant;
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
                        jinvariant.resize(16, '0');
                        if (tmp.size() > 1)
                            snprintf(jinvariant.data(), 17, "%08X%08X", tmp.data()[1], tmp.data()[0]);
                        else if (tmp.size() > 0)
                            snprintf(jinvariant.data() + 8, 9, "%08X", tmp.data()[0]);
                        dhash = tmp.to_string();
                        P.serialize(X, Y, Z, T);
                        EdD = ed_d;
                    }

                    std::unique_ptr<EdECMParams> params_edecm;
                    params_edecm.reset(new EdECMParams(net.task()->b1, net.task()->b2, maxSize, poly, polyThreads, polyMemModel));
                    logging.progress().add_stage(params_edecm->stage1_cost());
                    logging.progress().add_stage(params_edecm->stage2_cost());
                    net.task()->W = params_edecm->W;
                    net.task()->D = params_edecm->D;
                    net.task()->degree = 1 << params_edecm->PolyPower;

                    uint32_t fingerprint = File::unique_fingerprint(gwstate.fingerprint, jinvariant + "." + std::to_string(params_edecm->B1));
                    NetFile& file1 = files.emplace_back(net, "checkpoint.ed1", File::unique_fingerprint(fingerprint, std::to_string(params_edecm->W)));
                    NetFile& file12 = files.emplace_back(net, "checkpoint.ed12", fingerprint);
                    NetFile& file2 = files.emplace_back(net, "checkpoint.ed2", File::unique_fingerprint(fingerprint, std::to_string(params_edecm->B2)));
                    EdECMStage1::State* interstate = read_state<EdECMStage1::State>(&file12);
                    if (interstate == nullptr)
                    {
                        EdECMStage1 stage1((int)params_edecm->B1, params_edecm->W);
                        stage1.init(&input, &gwstate, &file1, &logging, &X, &Y, &Z, &T, &EdD);
                        stage1.run();
                        if (!stage1.success() && params_edecm->B2 > params_edecm->B1)
                        {
                            interstate = new EdECMStage1::State(std::move(*stage1.state()));
                            file12.write(*interstate);
                        }
                    }
                    logging.progress().next_stage();
                    if (interstate != nullptr)
                    {
                        std::unique_ptr<Stage2> stage2;
                        if (params_edecm->PolyPower == 0)
                            stage2.reset(new EdECMStage2(params_edecm->B1, params_edecm->B2, params_edecm->D, params_edecm->L, params_edecm->LN, logging));
                        else
                            stage2.reset(new EdECMStage2Poly(params_edecm->B1, params_edecm->B2, params_edecm->D, params_edecm->PolyPower, polyThreads, polyMemModel, polyCheck > 0));
                        dynamic_cast<IEdECMStage2*>(stage2.get())->init(&input, &gwstate, &file2, &logging, interstate->X(), interstate->Y(), interstate->Z(), interstate->T(), EdD);
                        stage2->run();
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
                        .Argument("result", net.task()->factor)
                        .Argument("dhash", dhash)
                        .Argument("time", std::to_string(logging.progress().time_total()))
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

