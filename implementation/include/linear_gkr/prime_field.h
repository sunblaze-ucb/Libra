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
		field_element operator + (const field_element &b) const;
		field_element operator * (const field_element &b) const;
		field_element operator / (const field_element &b) const;
		field_element operator - (const field_element &b) const;
		void set_value(const mpz_class &);
		std::string to_string(int);
	};
}
#endif