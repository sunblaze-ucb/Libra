#ifndef __prover
#define __prover

#include "traditional_gkr/circuit_bc_clogc.h"
#include "linear_gkr/prime_field.h"
#include "linear_gkr/polynomial.h"
#include <cstring>
#include <utility>
#include <chrono>
#include <iostream>
class prover
{
public:
	prime_field::field_element v_u, v_v;
	int total_uv;
	blocked_circuit C;
	prime_field::field_element*** circuit_value;

	int sumcheck_layer_id, length_g, length_u, length_v, block_binary_length;
	prime_field::field_element alpha, beta;
	const prime_field::field_element *r_0, *r_1, *r_b_0;
	prime_field::field_element *one_minus_r_0, *one_minus_r_1, *one_minus_r_b_0;

	quadratic_poly *add_linear, *multV_linear, *addV_linear, *mult_linear;
	prime_field::field_element *V_mult_add;
	prime_field::field_element *beta_u;
	prime_field::field_element *beta_g_r0_fhalf, *beta_g_r0_shalf, *beta_g_r1_fhalf, *beta_g_r1_shalf, *beta_u_fhalf, *beta_u_shalf;
	prime_field::field_element *beta_g_sum;
	prime_field::field_element *block_beta_value, *block_v_value;
	prime_field::field_element *V_mult_add_copy;
	prime_field::field_element beta_blk_val;
	double total_time;
	void init_array(int, int blk_bit_length);
	void get_circuit(const blocked_circuit &from_verifier);
	prime_field::field_element* evaluate();
	void proof_init();
	void sumcheck_init(int layer_id, int bit_length_blk, int bit_length_g, int bit_length_u, int bit_length_v, const prime_field::field_element &,
		const prime_field::field_element &, const prime_field::field_element*, const prime_field::field_element*,
		const prime_field::field_element*, prime_field::field_element*, prime_field::field_element*, prime_field::field_element*);
	void sumcheck_phase0_init();
	void sumcheck_phase1_init(prime_field::field_element);
	void sumcheck_phase2_init(prime_field::field_element);
	cubic_poly sumcheck_phase0_update(prime_field::field_element, int);
	quadratic_poly sumcheck_phase1_update(prime_field::field_element, int);
	quadratic_poly sumcheck_phase2_update(prime_field::field_element, int);
	prime_field::field_element V_res(const prime_field::field_element*, const prime_field::field_element*, const prime_field::field_element*, const prime_field::field_element*, const prime_field::field_element*, int, int, int);
	std::pair<prime_field::field_element, prime_field::field_element> sumcheck_finalize(prime_field::field_element);
	void delete_self();
	prover()
	{
	}
	~prover();
};

#endif
