#pragma once

#include <vector>
#include "arithmetic.h"

#ifndef FILE_APPID
#define FILE_APPID 1
#endif // !FILE_APPID

class Writer
{
public:
    Writer() { }
    Writer(std::vector<char>& buffer) : _buffer(buffer) { }
    Writer(std::vector<char>&& buffer) : _buffer(std::move(buffer)) { }

    void write(int32_t value);
    void write(uint32_t value);
    void write(uint64_t value);
    void write(double value);
    void write(const std::string& value);
    void write(const arithmetic::Giant& value);
    void write(const char* ptr, int count);

    std::vector<char>& buffer() { return _buffer; }
    std::vector<char> hash();
    std::string hash_str();

private:
    std::vector<char> _buffer;
};

class Reader
{
public:
    Reader(char format_version, char type, char version, char *data, int size, int pos) : _format_version(format_version), _type(type), _version(version), _data(data), _size(size), _pos(pos) { }

    bool read(int32_t& value);
    bool read(uint32_t& value);
    bool read(uint64_t& value);
    bool read(double& value);
    bool read(std::string& value);
    bool read(arithmetic::Giant& value);

    char type() { return _type; }
    char version() { return _version; }

private:
    char _format_version;
    char _type;
    char _version;
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
    virtual ~File() { }

    virtual Writer* get_writer();
    virtual Reader* get_reader();
    virtual void commit_writer(Writer& writer);
    virtual void clear();
    bool read(TaskState& state);
    void write(TaskState& state);
    
    std::string& filename() { return _filename; }
    bool hash = true;
    int appid = FILE_APPID;

protected:
    std::string _filename;
    uint32_t _fingerprint;
    std::vector<char> _buffer;
};
