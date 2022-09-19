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
	prime_field::field_element *beta_g_r0, *beta_g_r1, *beta_u, *beta_v;
	layered_circuit C;
	prover *p;
	void beta_init(int depth, prime_field::field_element alpha, prime_field::field_element beta,
	const prime_field::field_element* r_0, const prime_field::field_element* r_1, 
	const prime_field::field_element* r_u, const prime_field::field_element* r_v, 
	const prime_field::field_element* one_minus_r_0, const prime_field::field_element* one_minus_r_1, 
	const prime_field::field_element* one_minus_r_u, const prime_field::field_element* one_minus_r_v);
	void read_circuit(const char *);
	bool verify();
	void get_prover(prover*);
	void delete_self();
	prime_field::field_element mult(int);
	prime_field::field_element add(int);
	prime_field::field_element V_in(const prime_field::field_element*, const prime_field::field_element*, prime_field::field_element*, int, int);
	prime_field::field_element direct_relay(int depth, prime_field::field_element *r_g, prime_field::field_element *r_u);
	prime_field::field_element not_gate(int depth);
	prime_field::field_element minus_gate(int depth);
	prime_field::field_element xor_gate(int depth);
	prime_field::field_element NAAB_gate(int depth);
	prime_field::field_element sum_gate(int depth);
	prime_field::field_element relay_gate(int depth);
	void read_circuit_from_string(char*);
};

#endif
