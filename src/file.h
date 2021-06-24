#pragma once

#include <vector>
#include "arithmetic.h"

#ifndef FILE_APPID
#define FILE_APPID 1
#endif // !FILE_APPID

class Writer
{
public:
    Writer() { _data.reserve(_max_size); }

    void write(int32_t value);
    void write(uint32_t value);
    void write(double value);
    void write(const std::string& value);
    void write(const arithmetic::Giant& value);

    std::vector<char>& data() { return _data; }
    std::vector<char> hash();
    std::string hash_str();

private:
    std::vector<char> _data;
    static int _max_size;
};

class Reader
{
public:
    Reader(char format_version, char version, char *data, int size, int pos) : _format_version(format_version), _version(version), _data(data), _size(size), _pos(pos) { }
    Reader(char format_version, char version, std::vector<char>&& container, int pos) : _format_version(format_version), _version(version), _container(std::move(container)), _data(_container.data()), _size(_container.size()), _pos(pos) { }

    bool read(int32_t& value);
    bool read(uint32_t& value);
    bool read(double& value);
    bool read(std::string& value);
    bool read(arithmetic::Giant& value);

    char version() { return _version; }

private:
    char _format_version;
    char _version;
    std::vector<char> _container;
    char* _data;
    int _size;
    int _pos = 0;
};

class TaskState;

class File
{
public:
    static const uint32_t MAGIC_NUM = 0x9f2b3cd4;

public:
    File(const std::string& filename, uint32_t fingerprint) : _filename(filename), _fingerprint(fingerprint) { }

    virtual Reader* get_reader();
    virtual void commit_writer(Writer& writer);
    bool read(TaskState& state);
    void write(TaskState& state);
    void clear();
    
    bool hash = true;

private:
    std::string _filename;
    uint32_t _fingerprint;
};
