#include "linear_gkr/prime_field.h"

namespace prime_field
{
	mpz_class mod;
	bool initialized = false;
	void init(std::string s, int base)
	{
		assert(!initialized);
		initialized = true;
		mod.set_str(s, base);
	}
	void init(mpz_class m)
	{
		assert(!initialized);
		mod = m;
		initialized = true;
	}
	field_element::field_element(){}
	field_element::field_element(const mpz_class &t)
	{
		assert(initialized);
		value = t;
	}
	field_element::field_element(const int x)
	{
		assert(initialized);
		value = mpz_class(x);
	}
	field_element field_element::add_non_mod(const field_element &b)
	{
		assert(initialized);
		field_element ret;
		ret.value = value + b.value;
		return ret;
	}
	field_element field_element::operator + (const field_element &b) const
	{
		assert(initialized);
		field_element ret;
		ret.value = (b.value + value) % mod;
		return ret;
	}
	field_element field_element::operator * (const field_element &b) const
	{
		assert(initialized);
		field_element ret;
		ret.value = (b.value * value) % mod;
		return ret;
	}
	field_element field_element::operator / (const field_element &b) const
	{
		assert(initialized);
		field_element ret, inv;
		mpz_invert(inv.value.get_mpz_t(), b.value.get_mpz_t(), mod.get_mpz_t());
		ret.value = (value * inv.value) % mod;
		return ret;
	}
	field_element field_element::operator - (const field_element &b) const
	{
		assert(initialized);
		field_element ret;
		ret.value = (value + mod - b.value) % mod;
		return ret;
	}
	void field_element::set_value(const mpz_class &x)
	{
		value = x;
	}
	std::string field_element::to_string(int base = 10)	
	{
		return value.get_str(base);
	}
	bool field_element::operator != (const field_element &b) const
	{
		return value != b.value;
	}
}