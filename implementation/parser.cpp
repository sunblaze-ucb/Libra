#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <regex>
#include <vector>
#include <cassert>
using namespace std;
const int Maxn = 3000000;
const int Maxm = 3000000;

class Graph
{
public:
	int head[Maxn], Next[Maxm], No[Maxm], tot;
	void init(int n)
	{
		tot = 2;
		for(int i = 0; i <= n + 10; ++i)
			head[i] = 0;
	}
	void add(int x, int y)
	{
		Next[tot] = head[x];
		head[x] = tot;
		No[tot] = y;
		tot++;
	}
}G;

char source_gate0[1000], source_gate1[1000];
string source_line;
int target_gate, src0, src1, tgt;
int visit[Maxn];

vector<int> input_gates;

int mx_depth = 0;

void dfs(int x, int depth = 1)
{
	if(visit[x])
	{
		return;
	}
	if(mx_depth < depth)
		mx_depth = depth;
	visit[x] = depth;
	for(int i = G.head[x]; i; i = G.Next[i])
	{
		int y = G.No[i];
		dfs(y, depth + 1);
	}
}

int main()
{
	memset(visit, 0, sizeof visit);
	ifstream sha256in ("SHA256_64.pws");
	
	regex add_gate("P V[0-9]+ = V[0-9]+ \\+ V[0-9]+ E");
	regex mult_gate("P V[0-9]+ = V[0-9]+ \\* V[0-9]+ E");
	regex mult_constant_gate("P V[0-9]+ = [0-9]+ \\* V[0-9]+ E");
	regex xor_gate("P V[0-9]+ = V[0-9]+ XOR V[0-9]+ E");
	regex naab_gate("P V[0-9]+ = V[0-9]+ NAAB V[0-9]+ E");
	regex minus_gate("P V[0-9]+ = V[0-9]+ minus V[0-9]+ E");
	regex not_gate("P V[0-9]+ = V[0-9]+ NOT V[0-9]+ E");
	regex constant_assign_gate("P V[0-9]+ = [0-9]+ E");
	regex input_gate("P V[0-9]+ = I[0-9]+ E");
	regex output_gate("P O[0-9]+ = V[0-9]+ E");
	regex pass_gate("P V[0-9]+ = V[0-9]+ PASS V[0-9]+ E");
	regex sum_gate("P V[0-9]+ = V[0-9]+ SS V[0-9]+ E");
	
	smatch base_match;
	G.init(3000000 - 101);
	int tot_gates = 0;
	int deg2_gates = 0;
	while(getline(sha256in, source_line))
	{
		if(std::regex_match(source_line, base_match, add_gate))
		{
			sscanf(source_line.c_str(), "P V%d = V%d + V%d E", &tgt, &src0, &src1);
			G.add(src0, tgt);
			G.add(src1, tgt);
			tot_gates++;
		}
		else if(std::regex_match(source_line, base_match, mult_gate))
		{
			sscanf(source_line.c_str(), "P V%d = V%d * V%d E", &tgt, &src0, &src1);
			G.add(src0, tgt);
			G.add(src1, tgt);
			tot_gates++;
			deg2_gates++;
		}
		else if(std::regex_match(source_line, base_match, mult_constant_gate))
		{
			sscanf(source_line.c_str(), "P V%d = %d + V%d E", &tgt, &src0, &src1);
			G.add(src1, tgt);
			tot_gates++;
		}
		else if(std::regex_match(source_line, base_match, xor_gate))
		{
			sscanf(source_line.c_str(), "P V%d = V%d XOR V%d E", &tgt, &src0, &src1);
			G.add(src0, tgt);
			G.add(src1, tgt);
			tot_gates++;
			deg2_gates++;
		}
		else if(std::regex_match(source_line, base_match, naab_gate))
		{
			sscanf(source_line.c_str(), "P V%d = V%d NAAB V%d E", &tgt, &src0, &src1);
			G.add(src0, tgt);
			G.add(src1, tgt);
			tot_gates++;
			deg2_gates++;
		}
		else if(std::regex_match(source_line, base_match, minus_gate))
		{
			sscanf(source_line.c_str(), "P V%d = V%d minus V%d E", &tgt, &src0, &src1);
			G.add(src0, tgt);
			G.add(src1, tgt);
			tot_gates++;
		}
		else if(std::regex_match(source_line, base_match, not_gate))
		{
			sscanf(source_line.c_str(), "P V%d = V%d NOT V%d E", &tgt, &src0, &src1);
			G.add(src0, tgt);
			G.add(src1, tgt);
			tot_gates++;
		}
		else if(std::regex_match(source_line, base_match, constant_assign_gate))
		{
			tot_gates++;
		}
		else if(std::regex_match(source_line, base_match, input_gate))
		{
			sscanf(source_line.c_str(), "P V%d = I%d E", &tgt, &src0);
			src0 += 500000;
			G.add(src0, tgt);
			input_gates.push_back(src0);
			tot_gates++;
		}
		else if(std::regex_match(source_line, base_match, output_gate))
		{
			sscanf(source_line.c_str(), "P O%d = V%d E", &tgt, &src0);
			tgt += 1000000;
			G.add(src0, tgt);
			tot_gates++;
		}
		else if(std::regex_match(source_line, base_match, pass_gate))
		{
			sscanf(source_line.c_str(), "P V%d = V%d PASS V%d E", &tgt, &src0, &src1);
			tgt += 1000000;
			G.add(src0, tgt);
			tot_gates++;
		}
		else if(std::regex_match(source_line, base_match, sum_gate))
		{
			tot_gates++;
		}
		else
		{
			cout << source_line << endl;
			assert(false);
		}
	}
	for(int i = 0; i < input_gates.size(); ++i)
	{
		dfs(input_gates[i]);
	}
	cout << "max depth = " << mx_depth << endl;
	cout << "Total gates = " << tot_gates << endl;
	cout << "Deg2 gates = " << deg2_gates << endl;
	sha256in.close();
	return 0;
}
