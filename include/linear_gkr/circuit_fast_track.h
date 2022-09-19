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
	bool is_parallel;
	int block_size;
	int log_block_size;
	int repeat_num;
	int log_repeat_num;
	layer()
	{
		is_parallel = false;
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
