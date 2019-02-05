#include <algorithm>
#include <cstdio>
#include <ctime>
#include <vector>
#include <utility>
#include <map>
#include <cassert>

using namespace std;

pair<int, int> update_info[64][8];
vector<pair<int, int> > info_vec;
map<pair<int, int>, int> info_vec_addr;
FILE *fd;

map<pair<int, pair<int, int> >, int> var_to_aux_addr;

int bit_length(int x)
{
	int ret = 0;
	while(x)
		x = x / 2, ret++;
	return ret;
}

void construct_input(int n_leaf)
{
	//aux inputs
	//aux input W 64 words per instance
	//aux input abcdefgh 134 words per instance
	//total 198 words per instance with padding
	//n_leaf * 2 - 1 instances
	//total 198 * (n_leaf * 2 - 1) words in aux input
	//real input length = n_leaf * (512 / 32) words
	
	//last 256 bits final answer as aux
	
	//total inputs as bits 198 * 32 * (n_leaf * 2 - 1) = 12672 * n_leaf - 6336 bits
	
	//instances id from n_leaf to n_leaf * 2 - 1 are input instances
	//root instance has id 0
	
	int total_inputs = 0;
	int nbits = bit_length(12672 * n_leaf - 6336);
	int padding_num = (1 << nbits) - (12672 * n_leaf - 6336);
	fprintf(fd, "%d ", (1 << nbits));
	
	for(int i = 0; i < 2 * n_leaf - 1; ++i)
	{
		for(int j = 0; j < 64; ++j)
		{
			var_to_aux_addr[make_pair(i, make_pair('W', j))] = total_inputs;
			for(int k = 0; k < 32; ++k)
			{
				fprintf(fd, "3 %d %d 0 ", total_inputs, rand() % 2);
				total_inputs++;
			}
		}
		for(int j = 0; j < info_vec.size(); ++j)
		{
			var_to_aux_addr[make_pair(i, info_vec[j])] = total_inputs;
			for(int k = 0; k < 32; ++k)
			{
				fprintf(fd, "3 %d %d 0 ", total_inputs, rand() % 2);
				total_inputs++;
			}
		}
	}
	assert(total_inputs == 12672 * n_leaf - 6336);
	for(int i = 0; i < padding_num; ++i)
	{
		fprintf(fd, "3 %d 0 0", total_inputs++);
	}
}

int consistency_check_addr, bit_ops_addr;

void construct_bit_ops(int n_leaf)
{
	fprintf(fd, "%d ", (n_leaf - 1) * 512 + );
	//checking the input & output consistency
	int n_instances = n_leaf * 2 - 1;
	int tot_gates = 0;
	consistency_check_addr = 0;
	for(int i = 0; i < n_leaf - 1; ++i)
	{
		for(int j = 0; j < 8; ++j)
		{
			int input_instance_A;
			input_instance_A = i * 2;
			auto outA_var = update_info[63][j];
			int addrA = var_to_aux_addr[make_pair(input_instance_A, outA_var)];
			int addr_self = var_to_aux_addr[make_pair(i, make_pair('W', j))];
			for(int k = 0; k < 32; ++k)
			{
				fprintf(fd, "8 %d %d %d ", tot_gates++, addrA + k, addr_self + k);
			}
		}
		for(int j = 0; j < 8; ++j)
		{
			int input_instance_B;
			input_instance_B = i * 2 + 1;
			auto outB_var = update_info[63][j];
			int addrB = var_to_aux_addr[make_pair(input_instance_B, outB_var)];
			int addr_self = var_to_aux_addr[make_pair(i, make_pair('W', j + 8))];
			for(int k = 0; k < 32; ++k)
			{
				fprintf(fd, "8 %d %d %d ", tot_gates++, addrB + k, addr_self + k);
			}
		}
	}
	
	//bit_ops
	bit_ops_addr = tot_gates;
	for(int i = 0; i < n_instances; ++i)
	{
		for(int j = 16; j < 64; ++j)
		{
			auto w15_r7 = rotate_instr(i, j - 15, 7);
			auto w15_r18 = rotate_instr(i, j - 15, 18);
			auto w15_r3 = rotate_instr(i, j - 15, 3);
			
		}
	}
}

void preprocessing()
{
	update_info[0][0] = make_pair('a', 0);
	update_info[0][1] = make_pair('b', 0);
	update_info[0][2] = make_pair('c', 0);
	update_info[0][3] = make_pair('d', 0);
	update_info[0][4] = make_pair('e', 0);
	update_info[0][5] = make_pair('f', 0);
	update_info[0][6] = make_pair('g', 0);
	update_info[0][7] = make_pair('h', 0);
	
	for(int i = 1; i < 64; ++i)
	{
		update_info[i][0] = make_pair('a', i);
		update_info[i][1] = update_info[i - 1][0];
		update_info[i][2] = update_info[i - 1][1];
		update_info[i][3] = update_info[i - 1][2];
		update_info[i][4] = make_pair('e', i);
		update_info[i][5] = update_info[i - 1][4];
		update_info[i][6] = update_info[i - 1][5];
		update_info[i][7] = update_info[i - 1][6];
	}
	for(int i = 0; i < 64; ++i)
		for(int j = 0; j < 8; ++j)
			info_vec.push_back(update_info[i][j]);
	sort(info_vec.begin(), info_vec.end());
	info_vec.resize(unique(info_vec.begin(), info_vec.end()) - info_vec.begin());
	for(int i = 0; i < info_vec.size(); ++i)
		info_vec_addr[info_vec[i]] = i;
}



int main(int argc, char** argv)
{
	int n_leaf;
	fd = fopen(argv[1], "w");
	sscanf(argv[2], "%d", &n_leaf);
	
	preprocessing();
	
	fprintf(fd, "%d\n", 4);

	construct_input(n_leaf);
	construct_bit_ops(n_leaf);
	
	
	fclose(fd);
	return 0;
}
