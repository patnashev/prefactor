#pragma once

//#define FACTORING
//#define NETPF

#define PREFACTOR_VERSION "0.10.0"
#define VERSION_BUILD "27"

inline void print_banner()
{
    printf("Prefactor version " PREFACTOR_VERSION "." VERSION_BUILD ", GWnum library version " GWNUM_VERSION);
#ifdef GMP
    arithmetic::GMPArithmetic* gmp = dynamic_cast<arithmetic::GMPArithmetic*>(&arithmetic::GiantsArithmetic::default_arithmetic());
    if (gmp != nullptr)
        printf(", GMP library version %s", gmp->version().data());
#endif
    printf("\n");
}
