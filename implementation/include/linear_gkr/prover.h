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
	std::vector<int> current_sumcheck_gates;
	std::vector<std::unordered_map<int, prime_field::field_element > > circuit_value;

	int sumcheck_layer_id, length_g, length_u, length_v;
	std::vector<prime_field::field_element> randomness_from_verifier;
	prime_field::field_element alpha, beta;
	std::vector<prime_field::field_element> r_0, r_1; //previous randomness

	std::unordered_map<int, linear_poly> mult_array, addV_array, add_array;
	std::unordered_map<int, linear_poly> V_mult, V_add;
	clock_t total_time;

	void DFS(std::unordered_map<int, linear_poly> &, int, int, prime_field::field_element, 
		prime_field::field_element, const int);
	void DFS_add(std::unordered_map<int, linear_poly> &, int, int, const int);
	void get_circuit(const layered_circuit &from_verifier);
	std::vector<std::pair<int, prime_field::field_element> > evaluate();
	void proof_init();
	void sumcheck_init(int layer_id, int bit_length_g, int bit_length_u, int bit_length_v, const prime_field::field_element &,
		const prime_field::field_element &, const std::vector<prime_field::field_element>&, const std::vector<prime_field::field_element>&);
	void sumcheck_phase1_init();
	void sumcheck_phase2_init();
	quadratic_poly sumcheck_phase1_update(prime_field::field_element);
	quadratic_poly sumcheck_phase2_update(prime_field::field_element);
	prime_field::field_element V_0(const std::vector<prime_field::field_element> &r_0, 
								std::vector<std::pair<int, prime_field::field_element> > output);
};

#endif