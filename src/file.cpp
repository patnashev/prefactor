
#include <stdio.h>
#include <stdlib.h>
#include "gwnum.h"
#include "file.h"
#include "md5.c"
#include "inputnum.h"
#include "task.h"
#ifdef _WIN32
#include "windows.h"
#endif

int Writer::_max_size = 0;

void Writer::write(int32_t value)
{
    _data.insert(_data.end(), (char*)&value, 4 + (char*)&value);
    if (_max_size < _data.size())
        _max_size = (int)_data.size();
}

void Writer::write(uint32_t value)
{
    _data.insert(_data.end(), (char*)&value, 4 + (char*)&value);
    if (_max_size < _data.size())
        _max_size = (int)_data.size();
}

void Writer::write(double value)
{
    _data.insert(_data.end(), (char*)&value, (char*)(&value + 1));
    if (_max_size < _data.size())
        _max_size = (int)_data.size();
}

void Writer::write(const std::string& value)
{
    int32_t len = (int)value.size();
    _data.insert(_data.end(), (char*)&len, 4 + (char*)&len);
    _data.insert(_data.end(), (char*)value.data(), (char*)(value.data() + value.size()));
    if (_max_size < _data.size())
        _max_size = (int)_data.size();
}

void Writer::write(const arithmetic::Giant& value)
{
    int32_t len = value.size();
    if (value < 0)
        len *= -1;
    _data.insert(_data.end(), (char*)&len, 4 + (char*)&len);
    _data.insert(_data.end(), (char*)value.data(), (char*)(value.data() + value.size()));
    if (_max_size < _data.size())
        _max_size = (int)_data.size();
}

std::vector<char> Writer::hash()
{
    std::vector<char> digest(16);
    MD5_CTX context;
    MD5Init(&context);
    MD5Update(&context, (unsigned char *)_data.data(), (unsigned int)_data.size());
    MD5Final((unsigned char *)digest.data(), &context);
    return digest;
}

std::string Writer::hash_str()
{
    char output[33];
    md5_raw_input(output, (unsigned char *)_data.data(), (unsigned int)_data.size());
    return output;
}

bool Reader::read(int32_t& value)
{
    if (_size < _pos + 4)
        return false;
    value = *(int32_t*)(_data + _pos);
    _pos += 4;
    return true;
}

bool Reader::read(uint32_t& value)
{
    if (_size < _pos + 4)
        return false;
    value = *(uint32_t*)(_data + _pos);
    _pos += 4;
    return true;
}

bool Reader::read(double& value)
{
    if (_size < _pos + 4)
        return false;
    value = *(double*)(_data + _pos);
    _pos += sizeof(double);
    return true;
}

bool Reader::read(std::string& value)
{
    if (_size < _pos + 4)
        return false;
    int len = *(int32_t*)(_data + _pos);
    _pos += 4;
    if (_size < _pos + len)
        return false;
    value.insert(value.end(), _data + _pos, _data + _pos + len);
    _pos += len;
    return true;
}

bool Reader::read(arithmetic::Giant& value)
{
    if (_size < _pos + 4)
        return false;
    int len = *(int32_t*)(_data + _pos);
    _pos += 4;
    if (_size < _pos + abs(len)*4)
        return false;
    value.arithmetic().init((uint32_t*)(_data + _pos), abs(len), value);
    if (len < 0)
        value.arithmetic().neg(value, value);
    _pos += abs(len)*4;
    return true;
}

Reader* File::get_reader()
{
    Writer data;

    FILE* fd = fopen(_filename.data(), "rb");
    if (!fd)
        return nullptr;
    fseek(fd, 0L, SEEK_END);
    int filelen = ftell(fd);
    data.data().resize(filelen);
    fseek(fd, 0L, SEEK_SET);
    filelen = (int)fread(data.data().data(), 1, filelen, fd);
    fclose(fd);
    if (filelen != data.data().size())
        return nullptr;

    if (hash)
    {
        std::string md5_filename = _filename + ".md5";
        fd = fopen(md5_filename.data(), "rb");
        if (fd)
        {
            char saved_hash[33];
            fread(saved_hash, 1, 32, fd);
            fclose(fd);
            saved_hash[32] = 0;
            if (saved_hash != data.hash_str())
                return nullptr;
        }
    }

    if (data.data().size() < 8)
        return nullptr;
    if (*(uint32_t*)data.data().data() != MAGIC_NUM)
        return nullptr;
    if (data.data()[4] != FILE_APPID)
        return nullptr;

    return new Reader(data.data()[5], data.data()[6], std::move(data.data()), 8);
}

bool writeThrough(char *filename, const void *buffer, size_t count)
{
    bool ret = true;
#ifdef _WIN32
    HANDLE hFile = CreateFileA(filename, GENERIC_WRITE, FILE_SHARE_READ, NULL, CREATE_ALWAYS, FILE_FLAG_WRITE_THROUGH, NULL);
    if (hFile == INVALID_HANDLE_VALUE)
        return false;
    DWORD written;
    if (!WriteFile(hFile, buffer, (DWORD)count, &written, NULL) || written != count)
        ret = false;
    CloseHandle(hFile);
#else
    FILE *fd = fopen(filename, "wb");
    if (!fd)
        return false;
    if (fwrite(buffer, 1, count, fd) != count)
        ret = false;
    fclose(fd);
#endif
    return ret;
}

void File::commit_writer(Writer& writer)
{
    std::string new_filename = _filename + ".new";
    if (!writeThrough(new_filename.data(), writer.data().data(), writer.data().size()))
    {
        remove(new_filename.data());
        return;
    }
    remove(_filename.data());
    rename(new_filename.data(), _filename.data());

    if (hash)
    {
        new_filename = _filename + ".md5";
        std::string hash = writer.hash_str();
        writeThrough(new_filename.data(), hash.data(), 32);
    }
}

bool File::read(TaskState& state)
{
    std::unique_ptr<Reader> reader(get_reader());
    if (!reader)
        return false;
    uint32_t fingerprint;
    if (!reader->read(fingerprint) || fingerprint != _fingerprint)
        return false;
    return state.read(*reader);
}

void File::write(TaskState& state)
{
    Writer writer;
    writer.write(MAGIC_NUM);
    writer.write(FILE_APPID + (0 << 8) + (state.version() << 16));
    writer.write(_fingerprint);
    state.write(writer);
    commit_writer(writer);
}

void File::clear()
{
    remove(_filename.data());
    std::string md5_filename = _filename + ".md5";
    remove(md5_filename.data());
}
