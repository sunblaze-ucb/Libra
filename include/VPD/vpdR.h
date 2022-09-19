#ifndef __vpdR
#define __vpdR

#include "bn.h"
#include <vector>
#include "gmp.h"
#include "gmpxx.h"

namespace vpdR
{
    void KeyGen(int d);
    void environment_init();
    mpz_class commit(bn::Ec1& digest, bn::Ec1& digesta, std::vector<mpz_class>& input);
    void prove(std::vector<mpz_class> r, mpz_class& ans, std::vector<mpz_class>& input, std::vector<bn::Ec1>& witness, std::vector<bn::Ec1>& witnessa, mpz_class r_f);
    bool verify(std::vector<mpz_class> r, bn::Ec1 digest, mpz_class ans, std::vector<bn::Ec1>& witness, std::vector<bn::Ec1>& witnessa);
    bool check_commit(bn::Ec1 digest, bn::Ec1 digesta);
}

#endif