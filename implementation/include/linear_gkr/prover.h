#ifndef __prover
#define __prover

#include "linear_gkr/circuit.h"
#include <vector>
#include "linear_gkr/prime_field.h"
#include "linear_gkr/polynomial.h"
#include <ctime>

class prover
{
public:
	layered_circuit C;
	std::vector<std::unordered_map<int, prime_field::field_element > > circuit_value;

	int sumcheck_layer_id, length_g, length_u, length_v;
	std::vector<prime_field::field_element> randomness_from_verifier;
	prime_field::field_element alpha, beta;

	clock_t total_time;

	void get_circuit(const layered_circuit &from_verifier);
	std::vector<std::pair<int, prime_field::field_element> > evaluate();
	void proof_init();
	void sumcheck_init(int layer_id, int bit_length_g, int bit_length_u, int bit_length_v, const prime_field::field_element &,
		const prime_field::field_element &);
	void sumcheck_phase1_init();
	void sumcheck_phase2_init();
	prime_field::field_element V_0(const std::vector<prime_field::field_element> &r_0, 
								std::vector<std::pair<int, prime_field::field_element> > output);
};

#endif