
#include <cmath>
#include "primelist.h"

PrimeList::PrimeList(int max)
{
    int i, j, k;
    std::vector<char> bitmap(max/2, 0);
    k = 0;
    for (i = 1; i < max/2; i++)
        if (!bitmap[i])
        {
            k++;
            if (i > 16384)
                continue;
            for (j = (i*2 + 1)*(i*2 + 1)/2; j < bitmap.size(); j += i*2 + 1)
                bitmap[j] = 1;
        }
    primes.clear();
    primes.reserve(k + 1);
    primes.push_back(2);
    for (i = 1; i < max/2; i++)
        if (!bitmap[i])
            primes.push_back(i*2 + 1);
}

void PrimeList::sieve_range(int start, int end, std::vector<int>& list)
{
    int i, j;
    if (!(start & 1))
        start++;
    std::vector<char> bitmap((end - start)/2, 0);
    for (i = 1; i < primes.size() && primes[i]*primes[i] < end; i++)
    {
        j = (primes[i] - start%primes[i])%primes[i];
        for (j = ((j & 1) != 0 ? (j + primes[i]) : j)/2; j < bitmap.size(); j += primes[i])
            bitmap[j] = 1;
    }
    list.clear();
    list.reserve((int)((std::expint(log(end)) - std::expint(log(start)))*1.1) + 100);
    i = 0;
    if (start <= primes.back() && start < sqrt(end))
    {
        for (j = 0; primes[j] < start; j++);
        for (; j < primes.size() && primes[j]*primes[j] < end; j++)
            list.push_back(primes[j]);
    }
    for (j = 0; j < bitmap.size(); j++)
        if (!bitmap[j])
            list.push_back(start + j*2);
}

PrimeIterator PrimeList::begin()
{
    return PrimeIterator(*this);
}

PrimeIterator& PrimeIterator::operator++()
{
    _cur++;
    if (_cur < _list.size())
        return *this;
    if (_cur - _range_pos < _range.size())
        return *this;
    int start;
    if (_range_pos == 0)
    {
        _range_pos = _list.size();
        start = _list[_list.size() - 1] + 2;
    }
    else
    {
        _range_pos += _range.size();
        start = _range[_range.size() - 1] + 2;
    }
    _list.sieve_range(start, start + _list[_list.size() - 1], _range);
    return *this;
}

void PrimeIterator::operator++(int)
{
    operator++();
}

PrimeIterator& PrimeIterator::operator+=(int offset)
{
    _cur += offset;
    if (_cur < _list.size())
        return *this;
    while (_cur - _range_pos >= _range.size())
    {
        int start;
        if (_range_pos == 0)
        {
            _range_pos = _list.size();
            start = _list[_list.size() - 1] + 2;
        }
        else
        {
            _range_pos += _range.size();
            start = _range[_range.size() - 1] + 2;
        }
        _list.sieve_range(start, start + _list[_list.size() - 1], _range);
    }
    return *this;
}