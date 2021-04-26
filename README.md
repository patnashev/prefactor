# Prefactor
Prefactor is a command-line utility to find small factors of big numbers. It implements P-1, P+1 and EdECM methods of factoring using Gwnum library by George Woltman.

The main idea of all those methods is to build a sequence and check if it becomes periodic. It's possible to do it fast if the period is a multiple of several small primes, each of them not larger than a bound B. Such a period is called B-smooth. Its exact value depends on properies of the factor we're looking for and on properties of the method we use.

Probability of a number to be smooth depends on its size and known divisors. For example, if we know that *p-1* is divisible by Mersenne exponent, that makes the rest of *p-1* smaller and therefore more probable to be smooth. So, overall success of a factoring method depends among other things on how big are known divisors of the periods it is searching for.

All the methods use two-stage approach. The first stage is looking for B1-smooth periods, the second stage is looking for a period which is a multiple of primes all but one smaller than B1, and the last one smaller than B2. Usually B2 is 10-100 times larger than B1. Optimal B1 and B2 can be calculated by Prefactor from the sieving depth.

## P-1
P-1 is exceptionally good for Mersennes and GFNs, because *p-1* is divisible by 2p and 2^(n+1) respectively. For all other numbers P-1 is the worst, because *p-1* is only guaranteed to be divisible by 2 in general case.

## P+1
P+1 name is a little misleading, because it not always targets *p+1*. It targets either *p-1* or *p+1* depending on properties of *p*. So, if it is run after P-1, it can find only a half of factors, because the other half is already found. But P+1 takes in a parameter, which means P+1 can be run several times to cover most of *p+1*. But that strategy is not optimal for sieving.

There are special parameters of P+1, which guarantee that *p±1* has a divisor. *2/7* guarantees that *p±1* is divisible by 6. It "chooses" either *p-1* or *p+1* depending on which is divisible by 6. The same way *6/5* makes *p±1* divisible by 4. Those two parameters produce more factors than others, and more factors than P-1. It'd almost make sense to run P+1(*2/7*) as a preferred sieving method for general numbers, if not for its stage1 to be 50% slower compared to P-1. Higher probability of a factor does not compensate slower stage1. P-1 is still better as a sieve, it speeds up search for primes more.

## EdECM
Elliptic curve method (ECM) builds a sequence with *p+1+m* period. The main advantage of the method is that *m* is different for each curve we choose. If both *p-1* and *p+1* have large primes among their factors, we'll never find them with P-1 and P+1 methods. But with ECM we can run factorization multiple times with different curves and wait until we get *m* such that *p+1+m* is smooth enough for us to find it.

Operations with [Edwards curves](https://en.wikipedia.org/wiki/Edwards_curve) are faster than with any other curve type, that's why they're preferred for ECM.

As for the probability of success, curves can have known divisors of their periods too. There are several curve construction methods, the two most interesting for factoring guarantee that *p+1+m* is divisible by 12 and 16. It greatly increases smoothness probability, but not enough to compensate slow curve operations and make EdECM better than P-1 for sieving.

# Command line

```
Usage: prefactor {-B1 10000 -B2 100000 | -S sievingDepth [-B1 10000] [-B2 100000]} [-minus1] [-plus1] [-edecm] options {"K*B^N+C" | file}
Options: [-M maxMemory] [-t Threads] [-P 2/7] [-curve {curve2x8 | curve12 | random | seed 123 | xy 17/19 17/33}]
```
