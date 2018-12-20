#ifndef __circuit
#define __circuit
#include <utility>
#include "linear_gkr/prime_field.h"
#include <unordered_map>
#include <vector>

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
	gate *gates;
	int bit_length;
	std::unordered_map<int, std::vector<std::pair<int, std::pair<int, int> > > > u_gates;
	std::unordered_map<int, std::vector<std::pair<int, std::pair<int, int> > > > v_gates;
	layer()
	{
		gates = NULL;
		bit_length = 0;
	}
	~layer()
	{
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
		
	}
};

#endif
