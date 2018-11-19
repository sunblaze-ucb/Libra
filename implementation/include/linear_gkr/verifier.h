#pragma once
#ifndef __verifier
#define __verifier

#include "linear_gkr/circuit.h"
#include "linear_gkr/prover.h"
#include "linear_gkr/polynomial.h"

class verifier
{
public:
	layered_circuit C;
	prover *p;
	void read_circuit();
	bool verify();
	void get_prover(prover*);
};

#endif