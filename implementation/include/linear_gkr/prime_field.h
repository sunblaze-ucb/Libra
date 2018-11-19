#ifndef __prime_field
#define __prime_field

#include <gmp>
#include <gmpxx>

namespace prime_field
{
	gmp_class mod;
	bool initialized = false;
	void init(gmp_class m)
	{
		assert(!initialized);
		mod = m;
		initialized = true;
	}
	/*
	This defines a prime field
	*/
	class field_element
	{
	private:
		gmp_class value;
	public:
		field_element()
		{

		}
		~field_element()
		{

		}
		field_element operator + (const field_element &b) const
		{
			assert(initialized);
			field_element ret;
			ret.value = (b.value + value) % mod;
			return ret;
		}
		field_element operator * (const field_element &b) const
		{
			assert(initialized);
			field_element ret;
			ret.value = (b.value * value) % mod;
			return ret;
		}
		field_element operator / (const field_element &b) const
		{
			assert(initialized);
			field_element ret, inv;
			mpz_invert(inv.value,get_mpz_t(), b.value,get_mpz_t(), mod,get_mpz_t());
			ret.value = (value * inv);
			return ret;
		}
		field_element operator - (const field_element &b) const
		{
			assert(initialized);
			field_element ret;
			ret.value = (value + mod - b.value) % mod;
			return ret;
		}
		string to_string(int base = 10)
		{
			return value.get_str(base);
		}
	};
}
#endif