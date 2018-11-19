#include "linear_gkr/prover.h"

void prover::get_circuit(const layered_circuit &from_verifier)
{
	C = from_verifier;
}

std::vector<std::pair<int, prime_field::field_element> > prover::evaluate()
{
	circuit_value.clear();
	circuit_value.push_back(std::unordered_map<int, prime_field::field_element >());
	for(int i = 0; i < C.circuit[0].gate_id.size(); ++i)
	{
		int g, u, v, ty;
		g = C.circuit[0].gate_id[i];
		std::pair<int, std::pair<int, int> > info = C.circuit[0].gates[g];
		ty = info.first;
		u = info.second.first;
		v = info.second.second;
		if(ty == 0) //addition
		{
			circuit_value[0][g] = 
		}
		else if(ty == 1) //mult
		{

		}
		else if(ty == 2) //dummy
		{

		}
		else
		{
			assert(false);
		}
	}
}