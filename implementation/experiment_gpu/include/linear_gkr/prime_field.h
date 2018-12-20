#pragma once
#ifndef __prime_field
#define __prime_field

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/random.hpp>
#include <cassert>
#include <string>

using namespace boost::multiprecision;
using namespace boost::random;

namespace prime_field
{
	extern int512_t mod;
	extern bool initialized;
	extern independent_bits_engine<mt19937, 256, cpp_int> gen;
	void init(std::string, int);
	void init_random();
	/*
	This defines a prime field
	*/
	class field_element
	{
	private:
	public:
		int512_t value;

		field_element();
		field_element(const int);
		field_element operator + (const field_element &b) const;
		field_element operator * (const field_element &b) const;
		field_element operator / (const field_element &b) const;
		field_element operator - (const field_element &b) const;
		field_element mul_non_mod(const field_element &b) const;
		bool operator != (const field_element &b) const;
	};
	field_element random();
}
#endif