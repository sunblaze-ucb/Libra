#pragma once
#ifndef __verifier
#define __verifier

#include "linear_gkr/circuit_fast_track.h"
#include "linear_gkr/prover_fast_track.h"
#include "linear_gkr/polynomial.h"
#include <utility>

class verifier
{
public:
	layered_circuit C;
	prover *p;
	void read_circuit(const char *);
	bool verify();
	void get_prover(prover*);
	prime_field::field_element mult(int, const prime_field::field_element*, const prime_field::field_element*, const prime_field::field_element*, const prime_field::field_element*, const prime_field::field_element*, const prime_field::field_element*);
	prime_field::field_element add(int, const prime_field::field_element*, const prime_field::field_element*, const prime_field::field_element*, const prime_field::field_element*, const prime_field::field_element*, const prime_field::field_element*);
	prime_field::field_element V_in(const prime_field::field_element*, const prime_field::field_element*, prime_field::field_element*, int, int);
};

#endif