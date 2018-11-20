#include "linear_gkr/prover.h"

void prover::get_circuit(const layered_circuit &from_verifier)
{
	C = from_verifier;
}

prime_field::field_element prover::V_0(const std::vector<prime_field::field_element> &r_0, 
								std::vector<std::pair<int, prime_field::field_element> > output)
{
	std::vector<std::pair<int, prime_field::field_element> > tmp;
	for(int i = (int)r_0.size() - 1; i >= 0; --i)
	{
		tmp.clear();
		int last_gate;
		int cnt = 0;
		for(int j = 0; j < output.size(); ++j)
		{
			prime_field::field_element m = r_0[i];
			if((output[j].first & 1) == 0)
			{
				m = (prime_field::field_element(mpz_class(1)) - m);
			}
			if(j == 0)
			{
				tmp.push_back(std::make_pair(output[j].first >> 1, output[j].second * m));
				last_gate = output[j].first >> 1;
				cnt++;
			}
			else
			{
				if((output[j].first >> 1) == last_gate)
				{
					tmp[cnt - 1] = std::make_pair(last_gate, tmp[cnt - 1].second + output[j].second * m);
				}
				else
				{
					last_gate = output[j].first >> 1;
					tmp.push_back(std::make_pair(output[j].first >> 1, output[j].second * m));
					cnt++;
				}
			}
		}
		output = tmp;
	}
	assert(output.size() == 1);
	return output[0].second;
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
		assert(ty == 3);
		circuit_value[0][g] = mpz_class(v);
	}
	std::vector<std::pair<int, prime_field::field_element> > ret;
	for(int i = 1; i < C.circuit.size(); ++i)
	{
		circuit_value.push_back(std::unordered_map<int, prime_field::field_element>());
		for(int j = 0; j < C.circuit[i].gate_id.size(); ++j)
		{
			int g, u, v, ty;
			g = C.circuit[i].gate_id[j];
			std::pair<int, std::pair<int, int> > info = C.circuit[i].gates[g];
			ty = info.first;
			u = info.second.first;
			v = info.second.second;
			if(ty == 0)
			{
				circuit_value[i][g] = circuit_value[i - 1][u] + circuit_value[i - 1][v];
			}
			else if(ty == 1)
			{
				circuit_value[i][g] = circuit_value[i - 1][u] * circuit_value[i - 1][v];
			}
			else if(ty == 2)
			{
				circuit_value[i][g] = mpz_class(0);
			}
			else if(ty == 3)
			{
				circuit_value[i][g] = mpz_class(u);
			}
			else
			{
				assert(false);
			}
			if(i + 1 == C.circuit.size())
			{
				ret.push_back(std::make_pair(g, circuit_value[i][g]));
			}
		}
	}
	return ret;
}

void prover::proof_init()
{
	//todo
}