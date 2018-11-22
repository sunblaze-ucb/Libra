#pragma once
#ifndef __verifier
#define __verifier

#include "linear_gkr/circuit.h"
#include "linear_gkr/prover.h"
#include "linear_gkr/polynomial.h"

class verifier
{
public:
	layered_circuit C;
	prover *p;
	void read_circuit(const char *);
	bool verify();
	void get_prover(prover*);
	prime_field::field_element mult(int, std::vector<prime_field::field_element>, std::vector<prime_field::field_element>, std::vector<prime_field::field_element>);
	prime_field::field_element add(int, std::vector<prime_field::field_element>, std::vector<prime_field::field_element>, std::vector<prime_field::field_element>);
	prime_field::field_element V_in(const std::vector<prime_field::field_element> &, std::vector<std::pair<int, prime_field::field_element> >);
};

#endif