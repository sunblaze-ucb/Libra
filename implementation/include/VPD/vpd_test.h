#ifndef __vpd_test
#define __vpd_test

#include <cstdlib>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <gmp.h>
#include <gmpxx.h>
#include <math.h>
#include <string>

#include "test_point.hpp"
#include "bn.h"

using namespace std;
using namespace bn;

namespace vpd_test
{
    void KeyGen(int d);
    void environment_init();
    mpz_class commit(Ec1& digest, Ec1& digesta, vector<mpz_class>& input);
    bool check_commit(Ec1 digest, Ec1 digesta);
    void prove(vector<mpz_class> r, mpz_class& ans, vector<mpz_class>& input, vector<Ec1>& witness, vector<Ec1>& witnessa, mpz_class r_f);
    bool verify(vector<mpz_class> r, Ec1 digest, mpz_class& ans, vector<Ec1>& witness, vector<Ec1>& witnessa);
}

#endif