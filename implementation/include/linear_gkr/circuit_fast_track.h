#ifndef __circuit
#define __circuit
#include <utility>
#include "linear_gkr/prime_field.h"

class gate
{
public:
	int ty, u, v;
	gate(){}
	gate(int t, int U, int V)
	{
		ty = t, u = U, v = V;
	}
};

class layer
{
public:
	prime_field::field_element (*add)(std::vector<prime_field::field_element>);
	prime_field::field_element (*mult)(std::vector<prime_field::field_element>);
	gate *gates;
	int bit_length;
	layer()
	{
		gates = NULL;
		bit_length = 0;
	}
	~layer()
	{
		delete[] gates;
	}
};

class layered_circuit
{
public:
	layer *circuit;
	int total_depth;
	layered_circuit()
	{
		circuit = NULL;
		total_depth = 0;
	}
	~layered_circuit()
	{
		for(int i = 0; i < total_depth; ++i)
		{
			delete[] circuit[i].gates;
		}
		delete[] circuit;
	}
};

#endif