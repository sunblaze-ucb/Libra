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
	void get_circuit(const layered_circuit &from_verifier);
	std::vector<prime_field::field_element> evaluate(std::vector<prime_field::field_element> &input);
	void proof_init();
	polynomial sumcheck_step();
	polynomial sumcheck_final();
};

#endif