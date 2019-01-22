#ifndef __input_vpd
#define __input_vpd

#include <cstdlib>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <gmp.h>
#include <gmpxx.h>
#include <math.h>
#include <string>

#include "test_point.hpp"
#include "bn.h"

using namespace std;
using namespace bn;

namespace input_vpd
{
    mpz_class commit(Ec1& digest, Ec1& digesta, vector<mpz_class>& input);
    void KeyGen(int d);
    std::vector<mpz_class> pre_input(std::vector<mpz_class>& input);
    void prove(vector<mpz_class> r, mpz_class& ans, vector<mpz_class>& input, vector<Ec1>& witness, vector<Ec1>& witnessa, mpz_class r_f);
    void environment_init();
    bool verify(vector<mpz_class> r, Ec1 digest, mpz_class& ans, vector<Ec1>& witness, vector<Ec1>& witnessa);
}

#endif