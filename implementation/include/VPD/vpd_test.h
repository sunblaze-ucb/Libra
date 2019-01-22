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
}

#endif