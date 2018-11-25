#include "linear_gkr/prime_field.h"

namespace prime_field
{
	int512_t mod;
	bool initialized = false;
	independent_bits_engine<mt19937, 256, cpp_int> gen;
	void init(std::string s, int base)
	{
		assert(!initialized);
		assert(base == 10);
		initialized = true;
		mod = int512_t(s);
	}
	void init_random()
	{
	}
	field_element::field_element()
	{
		value = 0;
	}
	field_element::field_element(const int x)
	{
		assert(initialized);
		value = x;
	}
	field_element field_element::operator + (const field_element &b) const
	{
		assert(initialized);
		field_element ret;
		ret.value = (b.value + value);
		if(ret.value >= mod)
			ret.value = ret.value - mod;
		return ret;
	}
	field_element field_element::mul_non_mod(const field_element &b) const
	{
		assert(initialized);
		field_element ret;
		ret.value = (b.value * value);
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
		//todo
		assert(false);
		//return ret;
	}
	field_element field_element::operator - (const field_element &b) const
	{
		assert(initialized);
		field_element ret;
		if(value > b.value)
			ret.value = value - b.value;
		else
			ret.value = value + mod - b.value;
		return ret;
	}
	field_element random()
	{
		field_element ret;
		ret.value = int512_t(gen());
		return ret;
	}
	bool field_element::operator != (const field_element &b) const
	{
		return value != b.value;
	}
}