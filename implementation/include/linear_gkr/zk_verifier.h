#pragma once
#ifndef __zk_verifier
#define __zk_verifier

#include "linear_gkr/circuit_fast_track.h"
#include "linear_gkr/zk_prover.h"
#include "linear_gkr/polynomial.h"
#include "VPD/vpdR.h"
#include "VPD/vpd_test.h"
#include <utility>

class zk_verifier
{
public:
	//prime_field::field_element *beta_g_r0, *beta_g_r1, *beta_u, *beta_v;
	prime_field::field_element *beta_g_r0_first_half, *beta_g_r0_second_half;
	prime_field::field_element *beta_g_r1_first_half, *beta_g_r1_second_half;
	prime_field::field_element *beta_u_first_half, *beta_u_second_half;
	prime_field::field_element *beta_v_first_half, *beta_v_second_half;

	layered_circuit C;
	zk_prover *p;
	void beta_init(int depth, prime_field::field_element alpha, prime_field::field_element beta,
	const prime_field::field_element* r_0, const prime_field::field_element* r_1, 
	const prime_field::field_element* r_u, const prime_field::field_element* r_v, 
	const prime_field::field_element* one_minus_r_0, const prime_field::field_element* one_minus_r_1, 
	const prime_field::field_element* one_minus_r_u, const prime_field::field_element* one_minus_r_v);
	void read_circuit(const char *);
	bool verify();
	void get_prover(zk_prover*);
	void delete_self();
	prime_field::field_element mult(int);
	prime_field::field_element add(int);
	prime_field::field_element direct_relay(int depth, prime_field::field_element *r_g, prime_field::field_element *r_u);
	prime_field::field_element not_gate(int depth);
	prime_field::field_element minus_gate(int depth);
	prime_field::field_element xor_gate(int depth);
	prime_field::field_element NAAB_gate(int depth);
	prime_field::field_element sum_gate(int depth);
	prime_field::field_element V_in(const prime_field::field_element*, const prime_field::field_element*, prime_field::field_element*, int, int);
};

#endif
