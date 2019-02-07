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
//	prime_field::field_element (*add)(std::vector<prime_field::field_element>);
//	prime_field::field_element (*mult)(std::vector<prime_field::field_element>);
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

class blocked_circuit
{
public:
	layered_circuit* blocks;
	int total_blocks;
	int total_blocks_binary_length;
	blocked_circuit()
	{
		blocks = NULL;
		total_blocks = 0;
	}
};

#endif
