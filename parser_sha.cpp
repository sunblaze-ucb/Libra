#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <regex>
#include <map>
#include <vector>
#include <queue>
#include <cassert>
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

enum gate_types
{
    add = 0,
    mult = 1,
    dummy = 2,
    input = 3,
    not_gate_id = 6, 
    minus_gate_id = 7,
    xor_gate_id = 8,
    naab_gate_id = 9
};

class gate
{
public:
    gate_types ty;
    int id;
    long long u, v;
    bool u_modified, v_modified;
    int lvl;
    vector<int> outputs; //pointer to array
    gate()
    {
        lvl = -1;
    }
    bool operator < (const gate &b) const
    {
        return id < b.id;
    }
};

vector<gate> gates;
vector<gate> circuit_block_gates;
map<int, int> id_to_addr;
int repeat_num, rdl_output_cnt = 0;
int max_rdl_id = 0;

vector<int> rdl_output_pointer;

void read_rdl(ifstream &rdl_in)
{
    gate zero_gate;
    zero_gate.id = -1;
    zero_gate.ty = input;
    zero_gate.u = 0;
    zero_gate.v = 0;
    gates.push_back(zero_gate);
    while(getline(rdl_in, source_line))
	{
		if(std::regex_match(source_line, base_match, constant_assign_gate))
		{
            sscanf(source_line.c_str(), "P V%d = %lld E", &tgt, &src0);
            gate g;
            g.ty = input;
            g.id = tgt;
            g.u = src0;
            g.v = 0;
            gates.push_back(g);
            id_to_addr[g.id] = gates.size() - 1;
            if(tgt > max_rdl_id)
                max_rdl_id = tgt;
		}
		else if(std::regex_match(source_line, base_match, input_gate))
		{
			sscanf(source_line.c_str(), "P V%d = I%lld E", &tgt, &src0);
            gate g;
            g.ty = input;
            g.id = tgt;
            g.u = rand(); //random input, can be changed to more interesting input
            g.v = 0;
            gates.push_back(g);
            id_to_addr[g.id] = gates.size() - 1;
            if(tgt > max_rdl_id)
                max_rdl_id = tgt;
		}
		else if(std::regex_match(source_line, base_match, output_gate))
		{
            //doesn't need to do anything
			continue;
		}
		else if(std::regex_match(source_line, base_match, pass_gate))
		{
            rdl_output_cnt++;
			sscanf(source_line.c_str(), "P V%d = V%lld PASS V%lld E", &tgt, &src0, &src1);
            gate g;
            g.id = tgt;
            g.ty = add;
            g.u = src0;
            g.v = zero_gate.id;
            gates[id_to_addr[g.u]].outputs.push_back(tgt);
            gates[id_to_addr[g.v]].outputs.push_back(tgt);
            id_to_addr[g.id] = gates.size();
            rdl_output_pointer.push_back(gates.size());
            gates.push_back(g);
            if(tgt > max_rdl_id)
                max_rdl_id = tgt;
		}
		else
		{
			cout << source_line << endl;
			assert(false);
		}
	}
    max_rdl_id++;
    assert(rdl_output_cnt / repeat_num * repeat_num == rdl_output_cnt);
}

map<int, int> circuit_input;

void read_circuit(ifstream &circuit_in)
{
    while(getline(circuit_in, source_line))
    {
        if(std::regex_match(source_line, base_match, add_gate))
        {
            sscanf(source_line.c_str(), "P V%d = V%lld + V%lld E", &tgt, &src0, &src1);
            int u, v;
            for(int i = 0; i < repeat_num; ++i)
            {
                if(circuit_input.find(src0) != circuit_input.end())
                {
                    int rdl_gate_ptr = rdl_output_pointer[circuit_input[src0] * repeat_num + i];
                    u = gates[rdl_gate_ptr].id;
                }
                else
                {
                    u = src0 * repeat_num + i + max_rdl_id;
                }
                if(circuit_input.find(src1) != circuit_input.end())
                {
                    int rdl_gate_ptr = rdl_output_pointer[circuit_input[src1] * repeat_num + i];
                    v = gates[rdl_gate_ptr].id;
                }
                else
                {
                    v = src1 * repeat_num + i + max_rdl_id;
                }
                
                gate g;
                g.id = tgt * repeat_num + i + max_rdl_id;
                g.ty = add;
                g.u = u;
                g.v = v;
                id_to_addr[g.id] = gates.size();
                gates.push_back(g);
                assert(id_to_addr.find(u) != id_to_addr.end());
                assert(id_to_addr.find(v) != id_to_addr.end());
                gates[id_to_addr[u]].outputs.push_back(g.id);
                gates[id_to_addr[v]].outputs.push_back(g.id);
            }
            
        }
        else if(std::regex_match(source_line, base_match, mult_gate))
        {
            sscanf(source_line.c_str(), "P V%d = V%lld * V%lld E", &tgt, &src0, &src1);

            int u, v;
            for(int i = 0; i < repeat_num; ++i)
            {
                if(circuit_input.find(src0) != circuit_input.end())
                {
                    int rdl_gate_ptr = rdl_output_pointer[circuit_input[src0] * repeat_num + i];
                    u = gates[rdl_gate_ptr].id;
                }
                else
                {
                    u = src0 * repeat_num + i + max_rdl_id;
                }
                if(circuit_input.find(src1) != circuit_input.end())
                {
                    int rdl_gate_ptr = rdl_output_pointer[circuit_input[src1] * repeat_num + i];
                    v = gates[rdl_gate_ptr].id;
                }
                else
                {
                    v = src1 * repeat_num + i + max_rdl_id;
                }
                
                gate g;
                g.id = tgt * repeat_num + i + max_rdl_id;
                g.ty = mult;
                g.u = u;
                g.v = v;
                id_to_addr[g.id] = gates.size();
                gates.push_back(g);
                assert(id_to_addr.find(u) != id_to_addr.end());
                assert(id_to_addr.find(v) != id_to_addr.end());
                gates[id_to_addr[u]].outputs.push_back(g.id);
                gates[id_to_addr[v]].outputs.push_back(g.id);
            }
            
        }
        else if(std::regex_match(source_line, base_match, constant_assign_gate))
		{
            sscanf(source_line.c_str(), "P V%d = %lld E", &tgt, &src0);
            circuit_input[tgt] = src0;
		}
        else if(std::regex_match(source_line, base_match, input_gate))
        {
            sscanf(source_line.c_str(), "P V%d = I%lld E", &tgt, &src0);
            circuit_input[tgt] = src0;
        }
        else if(std::regex_match(source_line, base_match, output_gate))
        {
            sscanf(source_line.c_str(), "P O%d = V%lld E", &tgt, &src0);
            continue;
        }
        else if(std::regex_match(source_line, base_match, xor_gate))
        {
        	sscanf(source_line.c_str(), "P V%d = V%lld XOR V%lld E", &tgt, &src0, &src1);
        	int u, v;
            for(int i = 0; i < repeat_num; ++i)
            {
                if(circuit_input.find(src0) != circuit_input.end())
                {
                    int rdl_gate_ptr = rdl_output_pointer[circuit_input[src0] * repeat_num + i];
                    u = gates[rdl_gate_ptr].id;
                }
                else
                {
                    u = src0 * repeat_num + i + max_rdl_id;
                }
                if(circuit_input.find(src1) != circuit_input.end())
                {
                    int rdl_gate_ptr = rdl_output_pointer[circuit_input[src1] * repeat_num + i];
                    v = gates[rdl_gate_ptr].id;
                }
                else
                {
                    v = src1 * repeat_num + i + max_rdl_id;
                }
                
                gate g;
                g.id = tgt * repeat_num + i + max_rdl_id;
                g.ty = xor_gate_id;
                g.u = u;
                g.v = v;
                id_to_addr[g.id] = gates.size();
                gates.push_back(g);
                assert(id_to_addr.find(u) != id_to_addr.end());
                assert(id_to_addr.find(v) != id_to_addr.end());
                gates[id_to_addr[u]].outputs.push_back(g.id);
                gates[id_to_addr[v]].outputs.push_back(g.id);
            }
        }
        else if(std::regex_match(source_line, base_match, naab_gate))
        {
        	sscanf(source_line.c_str(), "P V%d = V%lld NAAB V%lld E", &tgt, &src0, &src1);
        	int u, v;
            for(int i = 0; i < repeat_num; ++i)
            {
                if(circuit_input.find(src0) != circuit_input.end())
                {
                    int rdl_gate_ptr = rdl_output_pointer[circuit_input[src0] * repeat_num + i];
                    u = gates[rdl_gate_ptr].id;
                }
                else
                {
                    u = src0 * repeat_num + i + max_rdl_id;
                }
                if(circuit_input.find(src1) != circuit_input.end())
                {
                    int rdl_gate_ptr = rdl_output_pointer[circuit_input[src1] * repeat_num + i];
                    v = gates[rdl_gate_ptr].id;
                }
                else
                {
                    v = src1 * repeat_num + i + max_rdl_id;
                }
                
                gate g;
                g.id = tgt * repeat_num + i + max_rdl_id;
                g.ty = naab_gate_id;
                g.u = u;
                g.v = v;
                id_to_addr[g.id] = gates.size();
                gates.push_back(g);
                assert(id_to_addr.find(u) != id_to_addr.end());
                assert(id_to_addr.find(v) != id_to_addr.end());
                gates[id_to_addr[u]].outputs.push_back(g.id);
                gates[id_to_addr[v]].outputs.push_back(g.id);
            }
        }
        else if(std::regex_match(source_line, base_match, minus_gate))
        {
        	sscanf(source_line.c_str(), "P V%d = V%lld minus V%lld E", &tgt, &src0, &src1);
        	int u, v;
            for(int i = 0; i < repeat_num; ++i)
            {
                if(circuit_input.find(src0) != circuit_input.end())
                {
                    int rdl_gate_ptr = rdl_output_pointer[circuit_input[src0] * repeat_num + i];
                    u = gates[rdl_gate_ptr].id;
                }
                else
                {
                    u = src0 * repeat_num + i + max_rdl_id;
                }
                if(circuit_input.find(src1) != circuit_input.end())
                {
                    int rdl_gate_ptr = rdl_output_pointer[circuit_input[src1] * repeat_num + i];
                    v = gates[rdl_gate_ptr].id;
                }
                else
                {
                    v = src1 * repeat_num + i + max_rdl_id;
                }
                
                gate g;
                g.id = tgt * repeat_num + i + max_rdl_id;
                g.ty = minus_gate_id;
                g.u = u;
                g.v = v;
                id_to_addr[g.id] = gates.size();
                gates.push_back(g);
                assert(id_to_addr.find(u) != id_to_addr.end());
                assert(id_to_addr.find(v) != id_to_addr.end());
                gates[id_to_addr[u]].outputs.push_back(g.id);
                gates[id_to_addr[v]].outputs.push_back(g.id);
            }
        }
        else if(std::regex_match(source_line, base_match, not_gate))
        {
        	sscanf(source_line.c_str(), "P V%d = V%lld NOT V%lld E", &tgt, &src0, &src1);
        	int u, v;
            for(int i = 0; i < repeat_num; ++i)
            {
                if(circuit_input.find(src0) != circuit_input.end())
                {
                    int rdl_gate_ptr = rdl_output_pointer[circuit_input[src0] * repeat_num + i];
                    u = gates[rdl_gate_ptr].id;
                }
                else
                {
                    u = src0 * repeat_num + i + max_rdl_id;
                }
                
                gate g;
                g.id = tgt * repeat_num + i + max_rdl_id;
                g.ty = not_gate_id;
                g.u = u;
                g.v = 0;
                id_to_addr[g.id] = gates.size();
                gates.push_back(g);
                assert(id_to_addr.find(u) != id_to_addr.end());
                gates[id_to_addr[u]].outputs.push_back(g.id);
            }
        }
        else
        {
            cout << source_line << endl;
            assert(false);
        }
    }
}

int max_lvl = -1;

map<int, pair<int, int> > id_to_layered_circuit_addr;

void top_sort()
{
    map<int, int> in_deg;
    for(int i = 0; i < gates.size(); ++i)
        in_deg[gates[i].id] = 0;
    for(int i = 0; i < gates.size(); ++i)
    {
        for(int j = 0; j < gates[i].outputs.size(); ++j)
            in_deg[gates[i].outputs[j]]++;
    }
    queue<int> q;
    for(int i = 0; i < gates.size(); ++i)
    {
        if(in_deg[gates[i].id] == 0)
            q.push(i), gates[i].lvl = 0;
    }
    while(!q.empty())
    {
        int cur = q.front();
        q.pop();
        for(int j = 0; j < gates[cur].outputs.size(); ++j)
        {
            int output_id = gates[cur].outputs[j];
            in_deg[output_id]--;
            int in_deg_v = in_deg[output_id];
            int output_addr = id_to_addr[output_id];
            gates[output_addr].lvl = max(gates[output_addr].lvl, gates[cur].lvl + 1);
            if(in_deg[output_id] == 0)
            {
                q.push(output_addr);
            }
        }
    }
    max_lvl = -1;
    for(int i = 0; i < gates.size(); ++i)
    {
        if(gates[i].lvl == -1)
        {
            printf("%d\n", in_deg[gates[i].id]);
        }
        assert(gates[i].lvl != -1);
        max_lvl = max(max_lvl, gates[i].lvl);
    }
}

int total_pad = 0, total_gates = 0;

void normalize_and_output()
{
    top_sort();
    vector<vector<gate> > layered_circuit;
    layered_circuit.resize(max_lvl + 1);
    //adding zero gate
    for(int i = 0; i < max_lvl + 1; ++i)
    {
        if(i == 0) //input layer zero gate
        {
            gate g;
            g.id = -10 - i;
            g.lvl = i;
            g.ty = input;
            g.u = 0;
            g.v = 0;
            id_to_addr[g.id] = gates.size();
            gates.push_back(g);
        }
        else
        {
            gate g;
            g.id = -10 - i;
            g.lvl = i;
            g.ty = add;
            g.u = -10 - (i - 1);
            g.v = -10 - (i - 1);
            id_to_addr[g.id] = gates.size();
            gates.push_back(g);
            gates[id_to_addr[g.u]].outputs.push_back(g.id);
            gates[id_to_addr[g.v]].outputs.push_back(g.id);
        }
    }
    int max_id = -1;
    for(int i = 0; i < gates.size(); ++i)
    {
        max_id = max(max_id, gates[i].id);
    }
    max_id = max_id + 1;
    //relay gate
    int cur_gates_sz = gates.size();
    for(int i = 0; i < cur_gates_sz; ++i)
    {
        for(int j = 0; j < gates[i].outputs.size(); ++j)
        {
            int output_id = gates[i].outputs[j];
            int output_addr = id_to_addr[output_id];
            int pre_gate_id = gates[i].id;
            bool need_relay = false;
            for(int lvl = gates[i].lvl + 1; lvl < gates[output_addr].lvl; ++lvl)
            {
                need_relay = true;
                gate g;
                g.id = max_id++;
                g.ty = add;
                g.u = pre_gate_id;
                g.v = -10 - (lvl - 1);
                g.lvl = lvl;
                if(lvl == gates[i].lvl + 1)
                {
                    gates[i].outputs[j] = g.id;
                }
                else
                    gates[id_to_addr[pre_gate_id]].outputs.push_back(g.id);
                gates[id_to_addr[g.v]].outputs.push_back(g.id);
                id_to_addr[g.id] = gates.size();
                gates.push_back(g);
                pre_gate_id = g.id;
            }
            if(gates[output_addr].u == gates[i].id && need_relay)
            {
                gates[output_addr].u = pre_gate_id;
                gates[id_to_addr[pre_gate_id]].outputs.push_back(gates[output_addr].id);
            }
            if(gates[output_addr].v == gates[i].id && need_relay)
            {
                gates[output_addr].v = pre_gate_id;
                gates[id_to_addr[pre_gate_id]].outputs.push_back(gates[output_addr].id);
            }
        }
    }
    for(int i = 0; i < gates.size(); ++i)
        layered_circuit[gates[i].lvl].push_back(gates[i]);

    
    for(int i = 0; i < layered_circuit.size(); ++i)
        sort(layered_circuit[i].begin(), layered_circuit[i].end());
    for(int i = 0; i < layered_circuit.size(); ++i)
    {
        for(int j = 0; j < layered_circuit[i].size(); ++j)
        {
            int id = layered_circuit[i][j].id;
            id_to_layered_circuit_addr[id] = make_pair(i, j);
            layered_circuit[i][j].u_modified = false;
            layered_circuit[i][j].v_modified = false;
        }
    }
    for(int i = 0; i < layered_circuit.size(); ++i)
    {
        for(int j = 0; j < layered_circuit[i].size(); ++j)
        {
            int new_id = j, old_id = layered_circuit[i][j].id;
            layered_circuit[i][j].id = new_id;
            for(int k = 0; k < layered_circuit[i][j].outputs.size(); ++k)
            {
                int output_id = layered_circuit[i][j].outputs[k];
                assert(id_to_layered_circuit_addr.find(output_id) != id_to_layered_circuit_addr.end());
                auto addr = id_to_layered_circuit_addr[output_id];
                bool modified = false;
                if(layered_circuit[addr.first][addr.second].u == old_id && !layered_circuit[addr.first][addr.second].u_modified)
                    layered_circuit[addr.first][addr.second].u = new_id, layered_circuit[addr.first][addr.second].u_modified = true, modified = true;
                if(layered_circuit[addr.first][addr.second].v == old_id && !layered_circuit[addr.first][addr.second].v_modified)
                    layered_circuit[addr.first][addr.second].v = new_id, layered_circuit[addr.first][addr.second].v_modified = true, modified = true;
                //assert(modified);
            }
        }
        int current_layer_sz = layered_circuit[i].size();
        if(__builtin_popcount(current_layer_sz) != 1) //padding
        {
            int target_size = 1;
            while(target_size < current_layer_sz)
                target_size *= 2;
            int padding_num = target_size - current_layer_sz;
            cout << current_layer_sz << " " << padding_num << endl;
            total_pad += padding_num;
            for(int j = 0; j < padding_num; ++j)
            {
                gate g;
                g.id = j + current_layer_sz;
                g.ty = dummy;
                g.u = 0;
                g.v = 0;
                g.u_modified = true;
                g.v_modified = true;
                layered_circuit[i].push_back(g);
            }
        }
    }
    FILE *out;
    out = fopen("sha_circuit.txt", "w");

    fprintf(out, "%d\n", (int)layered_circuit.size());
    for(int i = 0; i < layered_circuit.size(); ++i)
    {
        fprintf(out, "%d ", (int)layered_circuit[i].size());
        total_gates += layered_circuit[i].size();
        for(int j = 0; j < layered_circuit[i].size(); ++j)
        {
            fprintf(out, "%d %d %015lld %lld ", layered_circuit[i][j].ty, layered_circuit[i][j].id, layered_circuit[i][j].u, layered_circuit[i][j].v);
        }
        fprintf(out, "\n");
    }

    fclose(out);
}

int main(int argc, char* argv[])
{
	ifstream circuit_in (argv[1]);
    ifstream rdl_in(argv[2]);
    sscanf(argv[3], "%d", &repeat_num);
	
    read_rdl(rdl_in);
    read_circuit(circuit_in);
    normalize_and_output();
    cout << total_pad << " " << total_gates << endl;
    rdl_in.close();
	circuit_in.close();
	return 0;
}
