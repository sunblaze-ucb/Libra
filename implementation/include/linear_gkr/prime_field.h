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
	public:
		mpz_class value;
		field_element();
		field_element(const mpz_class&);
		field_element(const int);
		field_element add_non_mod(const field_element &b);
		field_element operator + (const field_element &b) const;
		field_element operator * (const field_element &b) const;
		field_element operator / (const field_element &b) const;
		field_element operator - (const field_element &b) const;
		bool operator != (const field_element &b) const;
		void set_value(const mpz_class &);
		std::string to_string(int);
	};
}
#endif