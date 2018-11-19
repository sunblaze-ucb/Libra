#pragma once
#ifndef __prime_field
#define __prime_field

#include <gmp.h>
#include <gmpxx.h>
#include <cassert>
#include <string>

namespace prime_field
{
	extern mpz_class mod;
	extern bool initialized;
	void init(std::string, int);
	void init(mpz_class);
	/*
	This defines a prime field
	*/
	class field_element
	{
	private:
		mpz_class value;
	public:
		inline field_element add_non_mod(const field_element &b);
		inline field_element operator + (const field_element &b) const;
		inline field_element operator * (const field_element &b) const;
		inline field_element operator / (const field_element &b) const;
		inline field_element operator - (const field_element &b) const;
		inline void set_value(const mpz_class &);
		inline std::string to_string(int);
	};
}
#endif