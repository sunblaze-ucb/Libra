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
    std::pair<mpz_class, mpz_class> commit(Ec1& digest, Ec1& digesta, Ec1& digest2, Ec1& digest2a, vector<mpz_class>& input, vector<mpz_class>& input2);
    void KeyGen(int d);
    void pre_input(std::vector<mpz_class>& input);
    void prove(vector<mpz_class> r, mpz_class& ans, vector<mpz_class>& input, vector<mpz_class> &input2, vector<Ec1>& witness, vector<Ec1>& witnessa, mpz_class r_f, mpz_class r_f2, mpz_class Z);
    void environment_init();
    bool verify(vector<mpz_class> r, Ec1 digest, Ec1 digest2, mpz_class Z, mpz_class& ans, vector<Ec1>& witness, vector<Ec1>& witnessa);
}

#endif