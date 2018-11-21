#ifndef __polynomial
#define __polynomial

#include <vector>
#include "gmp.h"
#include "gmpxx.h"
#include "linear_gkr/prime_field.h"


//ax^2 + bx + c
class quadratic_poly
{
public:
	prime_field::field_element a, b, c;
	quadratic_poly();
	quadratic_poly(const prime_field::field_element&, const prime_field::field_element&, const prime_field::field_element&);
	quadratic_poly operator + (const quadratic_poly &) const;
	prime_field::field_element eval(const prime_field::field_element &) const;
};

//ax + b
class linear_poly
{
public:
	prime_field::field_element a, b;
	linear_poly();
	linear_poly(const prime_field::field_element &, const prime_field::field_element &);
	linear_poly operator + (const linear_poly &) const;
	quadratic_poly operator * (const linear_poly &) const;
	prime_field::field_element eval(const prime_field::field_element &) const;
};


#endif