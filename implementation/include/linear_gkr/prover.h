#ifndef __prover
#define __prover

#include "linear_gkr/circuit.h"
#include <vector>
#include "linear_gkr/prime_field.h"
#include "linear_gkr/polynomial.h"

class prover
{
public:
	layered_circuit C;
	std::vector<std::unordered_map<int, prime_field::field_element > > circuit_value;
	void get_circuit(const layered_circuit &from_verifier);
	std::vector<std::pair<int, prime_field::field_element> > evaluate();
	void proof_init();
	void sumcheck_init();
	quadratic_poly sumcheck_step();
	prime_field::field_element V_0(const std::vector<prime_field::field_element> &r_0, 
								std::vector<std::pair<int, prime_field::field_element> > output);
};

#endif