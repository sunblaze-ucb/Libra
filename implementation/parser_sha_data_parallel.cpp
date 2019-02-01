#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <regex>
#include <map>
#include <vector>
#include <queue>
#include <cassert>
#include <algorithm>
using namespace std;

char source_gate0[1000], source_gate1[1000];
string source_line;
int target_gate, tgt;
long long src0, src1;
vector<int> input_gates;

regex add_gate("P V[0-9]+ = V[0-9]+ \\+ V[0-9]+ E");
regex mult_gate("P V[0-9]+ = V[0-9]+ \\* V[0-9]+ E");
regex constant_assign_gate("P V[0-9]+ = [\\-]*[0-9]+ E");
regex input_gate("P V[0-9]+ = I[0-9]+ E");
regex output_gate("P O[0-9]+ = V[0-9]+ E");
regex pass_gate("P V[0-9]+ = V[0-9]+ PASS V[0-9]+ E");
regex xor_gate("P V[0-9]+ = V[0-9]+ XOR V[0-9]+ E");
regex minus_gate("P V[0-9]+ = V[0-9]+ minus V[0-9]+ E");
regex naab_gate("P V[0-9]+ = V[0-9]+ NAAB V[0-9]+ E");
regex not_gate("P V[0-9]+ = V[0-9]+ NOT V[0-9]+ E");

smatch base_match;
int repeat_num;

enum gate_types
{
    add = 0,
    mult = 1,
    dummy = 2,
    input = 3,
    not_gate_id = 6, 
    minus_gate_id = 7,
    xor_gate_id = 8,
    naab_gate_id = 9,
    output_gate_id = 10
};

class gate
{
public:
	int ty, g;
	long long u, v;
	gate(){}
	gate(int TY, int G, int U, int V)
	{
		ty = TY;
		g = G;
		u = U;
		v = V;
	}
};

class layer
{
public:
	vector<gate> gates;
	int bit_len;
	bool is_parallel;
	int block_size;
	int log_block_size;
};

class layered_circuit
{
public:
	vector<layer> layers;
	int depth;
};

class DAG_gate
{
public:
	pair<int, int> input0, input1, id;
	int layered_id;
	int layered_lvl;
	bool is_output;
	gate_types ty;
	vector<pair<int, int> > outputs;
};

class DAG_circuit
{
public:
	map<pair<int, int>, DAG_gate> circuit;
};

layered_circuit rdl, sha256, sha256_combined;
DAG_circuit sha256_dag, sha256_dag_copy;

void DAG_to_layered()
{
	vector<int> layer_gate_count;
	vector<int> padding_num;
	map<pair<int, int>, int> in_deg;
	vector<pair<int, int> > sample_gate;
	for(auto x = sha256_dag.circuit.begin(); x != sha256_dag.circuit.end(); ++x)
	{
		auto &id = x -> first;
		auto &g = x -> second;
		for(auto y : g.outputs)
		{
			in_deg[y]++;
		}
		g.layered_lvl = -1;
	}
	queue<pair<int, int> > q;
	for(auto x = sha256_dag.circuit.begin(); x != sha256_dag.circuit.end(); ++x)
	{
		auto &id = x -> first;
		if(in_deg[id] == 0)
		{
			assert((x -> second).ty == input);
		}
		(x -> second).layered_lvl = 0;
		q.push(id);
	}
	int max_lvl = -1;
	while(!q.empty())
	{
		auto cur = q.front();
		q.pop();
		auto &g = sha256_dag.circuit[cur];
		max_lvl = max(max_lvl, g.layered_lvl);
		for(auto x : g.outputs)
		{
			in_deg[x]--;
			sha256_dag.circuit[x].layered_lvl = max(sha256_dag.circuit[x].layered_lvl, g.layered_lvl + 1);
			if(in_deg[x] == 0)
			{
				q.push(x);
			}
		}
	}
	sample_gate.resize(max_lvl + 1);
	for(auto x = sha256_dag.circuit.begin(); x != sha256_dag.circuit.end(); ++x)
	{
		auto &id = x -> first;
		if((x -> second).is_output)
			(x -> second).layered_lvl = max_lvl;
		sample_gate[(x -> second).layered_lvl] = x -> first;
	}
	//zero gate
	for(int i = 1; i <= max_lvl; ++i)
	{
		DAG_gate g;
		g.ty = minus;
		g.input0 = sample_gate[i - 1];
		g.input1 = sample_gate[i - 1];
		g.layered_lvl = i;
		g.id = make_pair((int)'A', 0);
		sha256_dag.circuit[g.id] = g;
	}
	//relay gates
	sha256_dag_copy = sha256;
	for(auto x = sha256_dag_copy.circuit.begin(); x != sha256_dag_copy.circuit.end(); ++x)
	{
		auto &id = x -> first;
		auto &g = x -> second;
		switch(g -> ty)
		{
			case add:
			{
				break;
			}
			case mult:
			{
				break;
			}
			case input:
			{
				break;
			}
			case not_gate_id:
			{
				break;
			}
			case minus_gate_id:
			{
				break;
			}
			case xor_gate_id:
			{
				break;
			}
			case naab_gate_id:
			{
				break;
			}
			case output_gate_id:
			{
				break;
			}
		}
	}
	layer_gate_count.resize(max_lvl + 1);
	padding_num.resize(max_lvl + 1);
	for(int i = 0; i <= max_lvl; ++i)
		layer_gate_count[i] = 0, padding_num[i] = 0;
	for(auto x = sha256_dag.circuit.begin(); x != sha256_dag.circuit.end(); ++x)
	{
		auto &id = x -> first;
		auto &g = x -> second;
		layer_gate_count[g.layered_lvl]++;
	}
	
	sha256.depth = max_lvl + 1;
	sha256.layers.resize(sha256.depth);
	for(int i = 0; i <= max_lvl; ++i)
	{
		int lgc = layer_gate_count[i];
		int full_layer_count = 1;
		int bit_length = 0;
		while(lgc)
		{
			full_layer_count *= 2;
			lgc /= 2;
			bit_length++;
		}
		padding_num[i] = full_layer_count - lgc;
		sha256.layers[i].bit_len = bit_length;
	}
	
	
}

void read_circuit(ifstream &circuit_in)
{
	while(getline(circuit_in, source_line))
    {
        if(std::regex_match(source_line, base_match, add_gate))
        {
            sscanf(source_line.c_str(), "P V%d = V%lld + V%lld E", &tgt, &src0, &src1);
            DAG_gate g;
            g.is_output = false;
            g.ty = add;
            g.id = make_pair((int)'V', (int)tgt);
            g.input0 = make_pair((int)'V', (int)src0);
            g.input1 = make_pair((int)'V', (int)src1);
            sha256_dag.circuit[make_pair((int)'V', tgt)] = g;
            sha256_dag.circuit[g.input0].outputs.push_back(g.id);
            sha256_dag.circuit[g.input1].outputs.push_back(g.id);
        }
        else if(std::regex_match(source_line, base_match, mult_gate))
        {
            sscanf(source_line.c_str(), "P V%d = V%lld * V%lld E", &tgt, &src0, &src1);
            DAG_gate g;
            g.is_output = false;
            g.ty = mult;
            g.id = make_pair((int)'V', (int)tgt);
            g.input0 = make_pair((int)'V', (int)src0);
            g.input1 = make_pair((int)'V', (int)src1);
            sha256_dag.circuit[make_pair((int)'V', tgt)] = g;
            sha256_dag.circuit[g.input0].outputs.push_back(g.id);
            sha256_dag.circuit[g.input1].outputs.push_back(g.id);
        }
        else if(std::regex_match(source_line, base_match, constant_assign_gate))
		{
			assert(false);
            sscanf(source_line.c_str(), "P V%d = %lld E", &tgt, &src0);
		}
        else if(std::regex_match(source_line, base_match, input_gate))
        {
            sscanf(source_line.c_str(), "P V%d = I%lld E", &tgt, &src0);
            DAG_gate g;
            g.is_output = false;
            g.ty = input;
            g.id = make_pair((int)'V', (int)tgt);
            g.input0 = make_pair((int)'I', (int)src0);
            g.input1 = make_pair((int)'N', 0);
            sha256_dag.circuit[make_pair((int)'V', tgt)] = g;
        }
        else if(std::regex_match(source_line, base_match, output_gate))
        {
            sscanf(source_line.c_str(), "P O%d = V%lld E", &tgt, &src0);
            sha256_dag.circuit[make_pair((int)'V', (int)src0)].is_output = true;
        }
        else if(std::regex_match(source_line, base_match, xor_gate))
        {
        	sscanf(source_line.c_str(), "P V%d = V%lld XOR V%lld E", &tgt, &src0, &src1);
            DAG_gate g;
            g.is_output = false;
            g.ty = xor_gate_id;
            g.id = make_pair((int)'V', (int)tgt);
            g.input0 = make_pair((int)'V', (int)src0);
            g.input1 = make_pair((int)'V', (int)src1);
            sha256_dag.circuit[make_pair((int)'V', tgt)] = g;
            sha256_dag.circuit[g.input0].outputs.push_back(g.id);
            sha256_dag.circuit[g.input1].outputs.push_back(g.id);
        }
        else if(std::regex_match(source_line, base_match, naab_gate))
        {
        	sscanf(source_line.c_str(), "P V%d = V%lld NAAB V%lld E", &tgt, &src0, &src1);
            DAG_gate g;
            g.is_output = false;
            g.ty = naab_gate_id;
            g.id = make_pair((int)'V', (int)tgt);
            g.input0 = make_pair((int)'V', (int)src0);
            g.input1 = make_pair((int)'V', (int)src1);
            sha256_dag.circuit[make_pair((int)'V', tgt)] = g;
            sha256_dag.circuit[g.input0].outputs.push_back(g.id);
            sha256_dag.circuit[g.input1].outputs.push_back(g.id);
        }
        else if(std::regex_match(source_line, base_match, minus_gate))
        {
        	sscanf(source_line.c_str(), "P V%d = V%lld minus V%lld E", &tgt, &src0, &src1);
            DAG_gate g;
            g.is_output = false;
            g.ty = minus_gate_id;
            g.id = make_pair((int)'V', (int)tgt);
            g.input0 = make_pair((int)'V', (int)src0);
            g.input1 = make_pair((int)'V', (int)src1);
            sha256_dag.circuit[make_pair((int)'V', tgt)] = g;
            sha256_dag.circuit[g.input0].outputs.push_back(g.id);
            sha256_dag.circuit[g.input1].outputs.push_back(g.id);
        }
        else if(std::regex_match(source_line, base_match, not_gate))
        {
        	sscanf(source_line.c_str(), "P V%d = V%lld NOT V%lld E", &tgt, &src0, &src1);
            DAG_gate g;
            g.is_output = false;
            g.ty = not_gate_id;
            g.id = make_pair((int)'V', (int)tgt);
            g.input0 = make_pair((int)'V', (int)src0);
            g.input1 = make_pair((int)'V', (int)src1);
            sha256_dag.circuit[make_pair((int)'V', tgt)] = g;
            sha256_dag.circuit[g.input0].outputs.push_back(g.id);
            sha256_dag.circuit[g.input1].outputs.push_back(g.id);
        }
        else
        {
            cout << source_line << endl;
            assert(false);
        }
    }
    
	DAG_to_layered();
}

int main(int argc, char* argv[])
{
	ifstream circuit_in (argv[1]);
    ifstream rdl_in(argv[2]);
    sscanf(argv[3], "%d", &repeat_num);
	
    read_circuit(circuit_in);
    
    rdl_in.close();
	circuit_in.close();
	return 0;
}
