#include "linear_gkr/verifier.h"
#include <string>
#include <utility>

void verifier::get_prover(prover *pp)
{
	p = pp;
}

void verifier::read_circuit(const char *path)	
{
	int d;
	static char str[300];
	FILE *circuit_in;
	circuit_in = fopen(path, "r");

	fscanf(circuit_in, "%d", &d);
	int n;
	fscanf(circuit_in, "%d", &d);

	C.input_gates.clear();
	C.input_gate_id.clear();
	for(int i = 0; i < n; ++i)
	{
		int g;
		fscanf(circuit_in, "%d%s", &d, str);
		mpz_class t;
		t.set_str(std::string(str), 10);
		C.input_gates[d] = t;
		C.input_gate_id.push_back(d);
	}
	C.circuit.clear();
	for(int i = 0; i < d - 1; ++i)
	{
		C.circuit.push_back(layer());
		fscanf(circuit_in, "%d", &n);
		for(int j = 0; j < n; ++j)
		{
			int ty, g, u, v;
			fscanf(circuit_in, "%d%d%d%d", &ty, &g, &u, &v);
			C.circuit[i].gates[g] = std::make_pair(ty, std::make_pair(u, v));
			C.circuit[i].gate_id.push_back(g);
		}
	}


	fclose(circuit_in);
	assert(false);
	//TODO
}