#include "linear_gkr/verifier.h"

void verifier::get_prover(prover *pp)
{
	p = pp;
}

void verifier::read_circuit(char *path)	
{
	int d;
	FILE *circuit_in;
	circuit_in = fopen(path, "r");

	fscanf(circuit_in, "%d", &d);

	for(int i = 0; i < d - 1; ++i)
	{

	}


	fclose(circuit_in);
	assert(false);
	//TODO
}