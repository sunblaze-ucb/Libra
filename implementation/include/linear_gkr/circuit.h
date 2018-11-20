#ifndef __circuit
#define __circuit
#include <vector>
#include <unordered_map>
#include <utility>
#include "linear_gkr/prime_field.h"

class layer
{
public:
	prime_field::field_element (*add)(std::vector<prime_field::field_element>);
	prime_field::field_element (*mult)(std::vector<prime_field::field_element>);
	std::unordered_map<int, std::pair<int, std::pair<int, int> > > gates;
	std::vector<int> gate_id;
	int bit_length;
};

class layered_circuit
{
public:
	std::vector<layer> circuit;
};

#endif