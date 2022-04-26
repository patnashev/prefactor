#pragma once

#include <vector>
#include <iterator>

int phi(int N);

class PrimeIterator;

class PrimeList
{
    friend class PrimeIterator;

public:
    PrimeList(int max);

    void sieve_range(int start, int end, std::vector<int>& list);

    size_t size() const { return primes.size(); }
    int operator[] (size_t pos) const { return primes[pos]; }

    PrimeIterator begin();

    static PrimeList& primes_16bit() { if (!_list65536) _list65536.reset(new PrimeList(65536)); return *_list65536; }

private:
    std::vector<int> primes;

    static std::unique_ptr<PrimeList> _list65536;
};

class PrimeIterator : public std::iterator<std::input_iterator_tag, int, int, const int*, int>
{
public:
    PrimeIterator(PrimeList& list) : _list(list) { }
    PrimeIterator(const PrimeIterator& it) : _list(it._list), _cur(it._cur), _range(it._range), _range_pos(it._range_pos) { }
    PrimeIterator& operator=(const PrimeIterator& it)
    {
        _list = it._list; _cur = it._cur; _range = it._range; _range_pos = it._range_pos;
        return *this;
    }

    static PrimeIterator get() { return PrimeIterator(PrimeList::primes_16bit()); }

    void sieve_range(uint64_t start, uint64_t end, std::vector<uint64_t>& list);

    PrimeIterator& operator++();
    void operator++(int);
    PrimeIterator& operator+=(int offset);
    bool operator==(PrimeIterator other) const { return _cur == other._cur; }
    bool operator!=(PrimeIterator other) const { return !(*this == other); }
    int operator*() const { return _cur < _list.size() ? _list[_cur] : _range[_cur - _range_pos]; }
    size_t pos() const { return _cur; }

private:
    PrimeList& _list;
    size_t _cur = 0;
    std::vector<int> _range;
    size_t _range_pos = 0;
};
