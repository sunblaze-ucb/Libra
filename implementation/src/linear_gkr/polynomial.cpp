#include "linear_gkr/polynomial.h"

quadratic_poly::quadratic_poly(){}
quadratic_poly::quadratic_poly(const prime_field::field_element& aa, const prime_field::field_element& bb, const prime_field::field_element& cc)
{
	a = aa;
	b = bb;
	c = cc;
}

quadratic_poly quadratic_poly::operator + (const quadratic_poly &x) const
{
	return quadratic_poly(a + x.a, b + x.b, c + x.c);
}

prime_field::field_element quadratic_poly::eval(const prime_field::field_element &x) const
{
	return ((a * x) + b) * x + c;
}

linear_poly::linear_poly(){}
linear_poly::linear_poly(const prime_field::field_element& aa, const prime_field::field_element& bb)
{
	a = aa;
	b = bb;
}

linear_poly linear_poly::operator + (const linear_poly &x) const
{
	return linear_poly(a + x.a, b + x.b);
}

quadratic_poly linear_poly::operator * (const linear_poly &x) const
{
	return quadratic_poly(a * x.a, a * x.b + b * x.a, b * x.b);
}

prime_field::field_element linear_poly::eval(const prime_field::field_element &x) const
{
	return a * x + b;
}