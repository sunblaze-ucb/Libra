#include "linear_gkr/zk_verifier.h"
#include <string>
#include <utility>
#include <gmp.h>
#include <vector>
#include <gmpxx.h>
#include <iostream>
#include "linear_gkr/random_generator.h"

void zk_verifier::get_prover(zk_prover *pp)
{
	p = pp;
}

void zk_verifier::read_circuit(const char *path, const char *meta_path)
{
	int d;
	FILE *circuit_in;
	FILE *meta_in;

	meta_in = fopen(meta_path, "r");
	circuit_in = fopen(path, "r");

	fscanf(circuit_in, "%d", &d);
	int n;
	C.circuit = new layer[d + 1];
	C.total_depth = d + 1;
	int max_bit_length = -1;
	for(int i = 1; i <= d; ++i)
	{
		fscanf(circuit_in, "%d", &n);
		if(i != 1)
		{
			if(n == 1)
				C.circuit[i].gates = new gate[2];
			else
				C.circuit[i].gates = new gate[n];
		}
		else
		{
			if(n == 1)
			{
				C.circuit[0].gates = new gate[2];
				C.circuit[1].gates = new gate[2];
			}
			else
			{
				C.circuit[0].gates = new gate[n];
				C.circuit[1].gates = new gate[n];
			}
		}
		
		int max_gate = -1;
		int previous_g = -1;
		for(int j = 0; j < n; ++j)
		{
			int ty, g;
			long long u, v;
			fscanf(circuit_in, "%d%d%lld%lld", &ty, &g, &u, &v);
			if(ty != 3)
			{
				if(ty == 5)
				{
					assert(u >= 0 && u < (1 << C.circuit[i - 1].bit_length));
					assert(v > u && v <= (1 << C.circuit[i - 1].bit_length));
				}
				else
				{
					if(!(u >= 0 && u < (1 << C.circuit[i - 1].bit_length)))
						cout << ty << " " << g << " " << u << " " << v << " " << (1 << C.circuit[i - 1].bit_length) << endl;
					assert(u >= 0 && u < (1 << C.circuit[i - 1].bit_length));
					if(!(v >= 0 && v < (1 << C.circuit[i - 1].bit_length)))
						cout << ty << " " << g << " " << u << " " << v << " " << (1 << C.circuit[i - 1].bit_length) << endl;
					assert(v >= 0 && v < (1 << C.circuit[i - 1].bit_length));
				}
			}
			if(ty == 6)
			{
				if(v != 0)
					fprintf(stderr, "WARNING, v!=0 for NOT gate.\n");
				v = 0;
			}
			if(ty == 10)
			{
				if(v != 0)
					fprintf(stderr, "WARNING, v!=0 for relay gate. %d\n", i);
				v = 0;
			}
			if(ty == 13)
			{
				assert(u == v);
			}
			if(g != previous_g + 1)
			{
				printf("Error, gates must be in sorted order, and full [0, 2^n - 1]. %d %d %d %d\n", i, j, g, previous_g);
				exit(0);
			}
			previous_g = g;
			if(i != 1)
				C.circuit[i].gates[g] = gate(ty, u, v);
			else
			{
				assert(ty == 2 || ty == 3);
				C.circuit[1].gates[g] = gate(4, g, 0);
				C.circuit[0].gates[g] = gate(ty, u, v);
			}
		}
		max_gate = previous_g;
		int cnt = 0;
		while(max_gate)
		{
			cnt++;
			max_gate >>= 1;
		}
		max_gate = 1;
		while(cnt)
		{
			max_gate <<= 1;
			cnt--;
		}
		int mx_gate = max_gate;
		while(mx_gate)
		{
			cnt++;
			mx_gate >>= 1;
		}
		if(n == 1)
		{
			//add a dummy gate to avoid ill-defined layer.
			if(i != 1)
			{
				C.circuit[i].gates[max_gate] = gate(2, 0, 0);
				C.circuit[i].bit_length = cnt;
			}
			else
			{
				C.circuit[0].gates[max_gate] = gate(2, 0, 0);
				C.circuit[0].bit_length = cnt;
				C.circuit[1].gates[max_gate] = gate(4, 1, 0);
				C.circuit[1].bit_length = cnt;
			}
		}
		else
		{
			C.circuit[i].bit_length = cnt - 1;
			if(i == 1)
				C.circuit[0].bit_length = cnt - 1;
		}
		fprintf(stderr, "layer %d, bit_length %d\n", i, C.circuit[i].bit_length);
		if(C.circuit[i].bit_length > max_bit_length)
			max_bit_length = C.circuit[i].bit_length;
	}
	C.circuit[0].is_parallel = false;
	for(int i = 1; i <= d; ++i)
	{
		int is_para;
		fscanf(meta_in, "%d", &is_para);
		fscanf(meta_in, "%d%d%d%d", &C.circuit[i].block_size, &C.circuit[i].repeat_num, &C.circuit[i].log_block_size, &C.circuit[i].log_repeat_num);
		if(is_para)
			assert(1 << C.circuit[i].log_repeat_num == C.circuit[i].repeat_num);
		if(is_para)
		{
			C.circuit[i].is_parallel = true;
		}
		else
		{
			C.circuit[i].is_parallel = false;
		}
		
	}
	p -> init_array(max_bit_length);

	int first_half_len = max_bit_length / 2, second_half_len = max_bit_length - first_half_len;
	beta_g_r0_first_half = new prime_field::field_element[(1 << first_half_len)];
	beta_g_r0_second_half = new prime_field::field_element[(1 << second_half_len)];
	beta_g_r1_first_half = new prime_field::field_element[(1 << first_half_len)];
	beta_g_r1_second_half = new prime_field::field_element[(1 << second_half_len)];
	beta_v_first_half = new prime_field::field_element[(1 << first_half_len)];
	beta_v_second_half = new prime_field::field_element[(1 << second_half_len)];
	beta_u_first_half = new prime_field::field_element[(1 << first_half_len)];
	beta_u_second_half = new prime_field::field_element[(1 << second_half_len)];

	beta_g_r0_block_first_half = new prime_field::field_element[(1 << first_half_len)];
	beta_g_r0_block_second_half = new prime_field::field_element[(1 << second_half_len)];
	beta_g_r1_block_first_half = new prime_field::field_element[(1 << first_half_len)];
	beta_g_r1_block_second_half = new prime_field::field_element[(1 << second_half_len)];
	beta_v_block_first_half = new prime_field::field_element[(1 << first_half_len)];
	beta_v_block_second_half = new prime_field::field_element[(1 << second_half_len)];
	beta_u_block_first_half = new prime_field::field_element[(1 << first_half_len)];
	beta_u_block_second_half = new prime_field::field_element[(1 << second_half_len)];
	fclose(circuit_in);
	fclose(meta_in);
}

vector<prime_field::field_element> zk_verifier::predicates(int depth, prime_field::field_element *r_0, prime_field::field_element *r_1, prime_field::field_element *r_u, prime_field::field_element *r_v, prime_field::field_element alpha, prime_field::field_element beta)
{
	vector<prime_field::field_element> ret_para;
	vector<prime_field::field_element> ret;
	const int gate_type_count = 14;
	ret.resize(gate_type_count);
	ret_para.resize(gate_type_count);
	for(int i = 0; i < gate_type_count; ++i)
	{
		ret[i] = prime_field::field_element(0);
		ret_para[i] = prime_field::field_element(0);
	}
	if(depth == 1)
	{
		return ret;
	}
	bool debug_mode = false;
	if(C.circuit[depth].is_parallel)
	{
		int first_half_g = C.circuit[depth].log_block_size / 2;
		int second_half_g = C.circuit[depth].log_block_size - first_half_g;
		int first_half_uv = C.circuit[depth - 1].log_block_size / 2;
		int second_half_uv = C.circuit[depth - 1].log_block_size - first_half_uv;
		int block_size = C.circuit[depth].block_size;
		vector<prime_field::field_element> one_block_alpha, one_block_beta;
		one_block_alpha.resize(gate_type_count);
		one_block_beta.resize(gate_type_count);
		for(int i = 0; i < gate_type_count; ++i)
		{
			one_block_alpha[i] = prime_field::field_element(0);
			one_block_beta[i] = prime_field::field_element(0);
		}
		assert((1 << C.circuit[depth].log_block_size) == C.circuit[depth].block_size);
		for(int i = 0; i < (1 << C.circuit[depth].log_block_size); ++i)
		{
			int g = i, u = C.circuit[depth].gates[i].u, v = C.circuit[depth].gates[i].v;
			g = g & ((1 << C.circuit[depth].log_block_size) - 1);
			u = u & ((1 << C.circuit[depth - 1].log_block_size) - 1);
			v = v & ((1 << C.circuit[depth - 1].log_block_size) - 1);
			switch(C.circuit[depth].gates[i].ty)
			{
				case 0:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					int u_first_half = u & ((1 << first_half_uv) - 1);
					int u_second_half = u >> first_half_uv;
					int v_first_half = v & ((1 << first_half_uv) - 1);
					int v_second_half = v >> first_half_uv;
					auto uv_value = (beta_u_block_first_half[u_first_half].value * beta_u_block_second_half[u_second_half].value % prime_field::mod) * (beta_v_block_first_half[v_first_half].value * beta_v_block_second_half[v_second_half].value % prime_field::mod) % prime_field::mod;
					one_block_alpha[0].value = one_block_alpha[0].value + (beta_g_r0_block_first_half[g_first_half].value * beta_g_r0_block_second_half[g_second_half].value) % prime_field::mod * uv_value;
					one_block_beta[0].value = one_block_beta[0].value + (beta_g_r1_block_first_half[g_first_half].value * beta_g_r1_block_second_half[g_second_half].value) % prime_field::mod * uv_value;
					one_block_alpha[0].value = one_block_alpha[0].value % prime_field::mod;
					one_block_beta[0].value = one_block_beta[0].value % prime_field::mod;
					break;
				}
				case 1:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					int u_first_half = u & ((1 << first_half_uv) - 1);
					int u_second_half = u >> first_half_uv;
					int v_first_half = v & ((1 << first_half_uv) - 1);
					int v_second_half = v >> first_half_uv;
					auto uv_value = (beta_u_block_first_half[u_first_half].value * beta_u_block_second_half[u_second_half].value % prime_field::mod) * (beta_v_block_first_half[v_first_half].value * beta_v_block_second_half[v_second_half].value % prime_field::mod) % prime_field::mod;
					one_block_alpha[1].value = one_block_alpha[1].value + (beta_g_r0_block_first_half[g_first_half].value * beta_g_r0_block_second_half[g_second_half].value) % prime_field::mod * uv_value;
					one_block_beta[1].value = one_block_beta[1].value + (beta_g_r1_block_first_half[g_first_half].value * beta_g_r1_block_second_half[g_second_half].value) % prime_field::mod * uv_value;
					one_block_alpha[1].value = one_block_alpha[1].value % prime_field::mod;
					one_block_beta[1].value = one_block_beta[1].value % prime_field::mod;
					break;
				}
				case 2:
				{
					break;
				}
				case 3:
				{
					break;
				}
				case 4:
				{
					break;
				}
				case 5:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					
					auto beta_g_val_alpha = beta_g_r0_block_first_half[g_first_half] * beta_g_r0_block_second_half[g_second_half];
					auto beta_g_val_beta = beta_g_r1_block_first_half[g_first_half] * beta_g_r1_block_second_half[g_second_half];
					auto beta_v_0 = beta_v_block_first_half[0] * beta_v_block_second_half[0];
					for(int j = u; j < v; ++j)
					{
						int u_first_half = j & ((1 << first_half_uv) - 1);
						int u_second_half = j >> first_half_uv;
						one_block_alpha[5] = one_block_alpha[5] + beta_g_val_alpha * beta_v_0 * (beta_u_block_first_half[u_first_half] * beta_u_block_second_half[u_second_half]);
						one_block_beta[5] = one_block_beta[5] + beta_g_val_beta * beta_v_0 * (beta_u_block_first_half[u_first_half] * beta_u_block_second_half[u_second_half]);
					}
					one_block_alpha[5].value = one_block_alpha[5].value % prime_field::mod;
					one_block_beta[5].value = one_block_beta[5].value % prime_field::mod;
					break;
				}
				case 12:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					
					auto beta_g_val_alpha = beta_g_r0_block_first_half[g_first_half] * beta_g_r0_block_second_half[g_second_half];
					auto beta_g_val_beta = beta_g_r1_block_first_half[g_first_half] * beta_g_r1_block_second_half[g_second_half];
					auto beta_v_0 = beta_v_block_first_half[0] * beta_v_block_second_half[0];
					for(int j = u; j <= v; ++j)
					{
						int u_first_half = j & ((1 << first_half_uv) - 1);
						int u_second_half = j >> first_half_uv;
						one_block_alpha[12] = one_block_alpha[12] + beta_g_val_alpha * beta_v_0 * (beta_u_block_first_half[u_first_half] * beta_u_block_second_half[u_second_half]);
						one_block_beta[12] = one_block_beta[12] + beta_g_val_beta * beta_v_0 * (beta_u_block_first_half[u_first_half] * beta_u_block_second_half[u_second_half]);

						beta_v_0 = beta_v_0 + beta_v_0;
						if(beta_v_0.value >= prime_field::mod)
							beta_v_0.value = beta_v_0.value - prime_field::mod;
					}
					one_block_alpha[12].value = one_block_alpha[12].value % prime_field::mod;
					one_block_beta[12].value = one_block_beta[12].value % prime_field::mod;
					break;
				}
				case 6:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					int u_first_half = u & ((1 << first_half_uv) - 1);
					int u_second_half = u >> first_half_uv;
					int v_first_half = v & ((1 << first_half_uv) - 1);
					int v_second_half = v >> first_half_uv;
					auto uv_value = (beta_u_block_first_half[u_first_half].value * beta_u_block_second_half[u_second_half].value % prime_field::mod) * (beta_v_block_first_half[v_first_half].value * beta_v_block_second_half[v_second_half].value % prime_field::mod) % prime_field::mod;
					one_block_alpha[6].value = one_block_alpha[6].value + (beta_g_r0_block_first_half[g_first_half].value * beta_g_r0_block_second_half[g_second_half].value) % prime_field::mod * uv_value;
					one_block_beta[6].value = one_block_beta[6].value + (beta_g_r1_block_first_half[g_first_half].value * beta_g_r1_block_second_half[g_second_half].value) % prime_field::mod * uv_value;
					one_block_alpha[6].value = one_block_alpha[6].value % prime_field::mod;
					one_block_beta[6].value = one_block_beta[6].value % prime_field::mod;
					break;
				}
				case 7:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					int u_first_half = u & ((1 << first_half_uv) - 1);
					int u_second_half = u >> first_half_uv;
					int v_first_half = v & ((1 << first_half_uv) - 1);
					int v_second_half = v >> first_half_uv;
					auto uv_value = (beta_u_block_first_half[u_first_half].value * beta_u_block_second_half[u_second_half].value % prime_field::mod) * (beta_v_block_first_half[v_first_half].value * beta_v_block_second_half[v_second_half].value % prime_field::mod) % prime_field::mod;
					one_block_alpha[7].value = one_block_alpha[7].value + (beta_g_r0_block_first_half[g_first_half].value * beta_g_r0_block_second_half[g_second_half].value) % prime_field::mod * uv_value;
					one_block_beta[7].value = one_block_beta[7].value + (beta_g_r1_block_first_half[g_first_half].value * beta_g_r1_block_second_half[g_second_half].value) % prime_field::mod * uv_value;
					one_block_alpha[7].value = one_block_alpha[7].value % prime_field::mod;
					one_block_beta[7].value = one_block_beta[7].value % prime_field::mod;
					break;
				}
				case 8:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					int u_first_half = u & ((1 << first_half_uv) - 1);
					int u_second_half = u >> first_half_uv;
					int v_first_half = v & ((1 << first_half_uv) - 1);
					int v_second_half = v >> first_half_uv;
					auto uv_value = (beta_u_block_first_half[u_first_half].value * beta_u_block_second_half[u_second_half].value % prime_field::mod) * (beta_v_block_first_half[v_first_half].value * beta_v_block_second_half[v_second_half].value % prime_field::mod) % prime_field::mod;
					one_block_alpha[8].value = one_block_alpha[8].value + (beta_g_r0_block_first_half[g_first_half].value * beta_g_r0_block_second_half[g_second_half].value) % prime_field::mod * uv_value;
					one_block_beta[8].value = one_block_beta[8].value + (beta_g_r1_block_first_half[g_first_half].value * beta_g_r1_block_second_half[g_second_half].value) % prime_field::mod * uv_value;
					one_block_alpha[8].value = one_block_alpha[8].value % prime_field::mod;
					one_block_beta[8].value = one_block_beta[8].value % prime_field::mod;
					break;
				}
				case 9:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					int u_first_half = u & ((1 << first_half_uv) - 1);
					int u_second_half = u >> first_half_uv;
					int v_first_half = v & ((1 << first_half_uv) - 1);
					int v_second_half = v >> first_half_uv;
					auto uv_value = (beta_u_block_first_half[u_first_half].value * beta_u_block_second_half[u_second_half].value % prime_field::mod) * (beta_v_block_first_half[v_first_half].value * beta_v_block_second_half[v_second_half].value % prime_field::mod) % prime_field::mod;
					one_block_alpha[9].value = one_block_alpha[9].value + (beta_g_r0_block_first_half[g_first_half].value * beta_g_r0_block_second_half[g_second_half].value) % prime_field::mod * uv_value;
					one_block_beta[9].value = one_block_beta[9].value + (beta_g_r1_block_first_half[g_first_half].value * beta_g_r1_block_second_half[g_second_half].value) % prime_field::mod * uv_value;
					one_block_alpha[9].value = one_block_alpha[9].value % prime_field::mod;
					one_block_beta[9].value = one_block_beta[9].value % prime_field::mod;
					break;
				}
				case 10:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					int u_first_half = u & ((1 << first_half_uv) - 1);
					int u_second_half = u >> first_half_uv;
					int v_first_half = v & ((1 << first_half_uv) - 1);
					int v_second_half = v >> first_half_uv;
					auto uv_value = (beta_u_block_first_half[u_first_half].value * beta_u_block_second_half[u_second_half].value % prime_field::mod) * (beta_v_block_first_half[v_first_half].value * beta_v_block_second_half[v_second_half].value % prime_field::mod) % prime_field::mod;
					one_block_alpha[10].value = one_block_alpha[10].value + (beta_g_r0_block_first_half[g_first_half].value * beta_g_r0_block_second_half[g_second_half].value) % prime_field::mod * uv_value;
					one_block_beta[10].value = one_block_beta[10].value + (beta_g_r1_block_first_half[g_first_half].value * beta_g_r1_block_second_half[g_second_half].value) % prime_field::mod * uv_value;
					one_block_alpha[10].value = one_block_alpha[10].value % prime_field::mod;
					one_block_beta[10].value = one_block_beta[10].value % prime_field::mod;
					break;
				}
				case 13:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					int u_first_half = u & ((1 << first_half_uv) - 1);
					int u_second_half = u >> first_half_uv;
					int v_first_half = v & ((1 << first_half_uv) - 1);
					int v_second_half = v >> first_half_uv;
					auto uv_value = (beta_u_block_first_half[u_first_half].value * beta_u_block_second_half[u_second_half].value % prime_field::mod) * (beta_v_block_first_half[v_first_half].value * beta_v_block_second_half[v_second_half].value % prime_field::mod) % prime_field::mod;
					one_block_alpha[13].value = one_block_alpha[13].value + (beta_g_r0_block_first_half[g_first_half].value * beta_g_r0_block_second_half[g_second_half].value) % prime_field::mod * uv_value;
					one_block_beta[13].value = one_block_beta[13].value + (beta_g_r1_block_first_half[g_first_half].value * beta_g_r1_block_second_half[g_second_half].value) % prime_field::mod * uv_value;
					one_block_alpha[13].value = one_block_alpha[13].value % prime_field::mod;
					one_block_beta[13].value = one_block_beta[13].value % prime_field::mod;
					break;
				}
			}
		}
		for(int i = 0; i < C.circuit[depth].repeat_num; ++i)
		{
			prime_field::field_element prefix_alpha, prefix_beta;
			prime_field::field_element prefix_alpha_v0, prefix_beta_v0;
			prefix_alpha = prime_field::field_element(1);
			prefix_beta = prime_field::field_element(1);
			prefix_alpha_v0 = prime_field::field_element(1);
			prefix_beta_v0 = prime_field::field_element(1);
			for(int j = 0; j < C.circuit[depth].log_repeat_num; ++j)
			{
				if((i >> j) & 1)
				{
					auto uv_value = r_u[j + C.circuit[depth - 1].log_block_size] * r_v[j + C.circuit[depth - 1].log_block_size];
					prefix_alpha = prefix_alpha * r_0[j + C.circuit[depth].log_block_size] * uv_value;
					prefix_beta = prefix_beta * r_1[j + C.circuit[depth].log_block_size] * uv_value;
				}
				else
				{
					auto uv_value = (prime_field::field_element(1) - r_u[j + C.circuit[depth - 1].log_block_size]) * (prime_field::field_element(1) - r_v[j + C.circuit[depth - 1].log_block_size]);
					prefix_alpha = prefix_alpha * (prime_field::field_element(1) - r_0[j + C.circuit[depth].log_block_size]) * uv_value;
					prefix_beta = prefix_beta * (prime_field::field_element(1) - r_1[j + C.circuit[depth].log_block_size]) * uv_value;
				}
			}
			for(int j = 0; j < C.circuit[depth].log_repeat_num; ++j)
			{
				if((i >> j) & 1)
				{
					auto uv_value = r_u[j + C.circuit[depth - 1].log_block_size] * (prime_field::field_element(1) - r_v[j + C.circuit[depth - 1].log_block_size]);
					prefix_alpha_v0 = prefix_alpha_v0 * r_0[j + C.circuit[depth].log_block_size] * uv_value;
					prefix_beta_v0 = prefix_beta_v0 * r_1[j + C.circuit[depth].log_block_size] * uv_value;
				}
				else
				{
					auto uv_value = (prime_field::field_element(1) - r_u[j + C.circuit[depth - 1].log_block_size]) * (prime_field::field_element(1) - r_v[j + C.circuit[depth - 1].log_block_size]);
					prefix_alpha_v0 = prefix_alpha_v0 * (prime_field::field_element(1) - r_0[j + C.circuit[depth].log_block_size]) * uv_value;
					prefix_beta_v0 = prefix_beta_v0 * (prime_field::field_element(1) - r_1[j + C.circuit[depth].log_block_size]) * uv_value;
				}
			}
			for(int j = 0; j < gate_type_count; ++j)
			{
				if(j == 6 || j == 10 || j == 5 || j == 12)
				{
					ret_para[j].value = ret_para[j].value + prefix_alpha_v0.value * one_block_alpha[j].value + prefix_beta_v0.value * one_block_beta[j].value;
					ret_para[j].value = ret_para[j].value % prime_field::mod;
				}
				else
				{
					ret_para[j].value = ret_para[j].value + prefix_alpha.value * one_block_alpha[j].value + prefix_beta.value * one_block_beta[j].value;
					ret_para[j].value = ret_para[j].value % prime_field::mod;
				}
			}
		}
		if(!debug_mode)
			ret = ret_para;
	}
	if(!C.circuit[depth].is_parallel || debug_mode)
	{
		int first_half_g = C.circuit[depth].bit_length / 2;
		int second_half_g = C.circuit[depth].bit_length - first_half_g;
		int first_half_uv = C.circuit[depth - 1].bit_length / 2;
		int second_half_uv = C.circuit[depth - 1].bit_length - first_half_uv;

		prime_field::field_element *tmp_u_val = new prime_field::field_element[1 << C.circuit[depth - 1].bit_length];
		prime_field::field_element zero_v;
		zero_v.value = (beta_v_first_half[0].value * beta_v_second_half[0].value % prime_field::mod);
		for(int i = 0; i < (1 << C.circuit[depth - 1].bit_length); ++i)
		{
			int u_first_half = i & ((1 << first_half_uv) - 1);
			int u_second_half = i >> first_half_uv;
			tmp_u_val[i].value = (beta_u_first_half[u_first_half].value * beta_u_second_half[u_second_half].value % prime_field::mod) * zero_v.value % prime_field::mod;
		}

		for(int i = 0; i < (1 << C.circuit[depth].bit_length); ++i)
		{
			int g = i, u = C.circuit[depth].gates[i].u, v = C.circuit[depth].gates[i].v;
			switch(C.circuit[depth].gates[i].ty)
			{
				case 0:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					int u_first_half = u & ((1 << first_half_uv) - 1);
					int u_second_half = u >> first_half_uv;
					int v_first_half = v & ((1 << first_half_uv) - 1);
					int v_second_half = v >> first_half_uv;
					ret[0].value = ret[0].value + (beta_g_r0_first_half[g_first_half].value * beta_g_r0_second_half[g_second_half].value + beta_g_r1_first_half[g_first_half].value * beta_g_r1_second_half[g_second_half].value) % prime_field::mod * 
								(beta_u_first_half[u_first_half].value * beta_u_second_half[u_second_half].value % prime_field::mod) % prime_field::mod * (beta_v_first_half[v_first_half].value * beta_v_second_half[v_second_half].value % prime_field::mod);
					ret[0].value = ret[0].value % prime_field::mod;
					break;
				}
				case 1:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					int u_first_half = u & ((1 << first_half_uv) - 1);
					int u_second_half = u >> first_half_uv;
					int v_first_half = v & ((1 << first_half_uv) - 1);
					int v_second_half = v >> first_half_uv;
					ret[1].value = ret[1].value + (beta_g_r0_first_half[g_first_half].value * beta_g_r0_second_half[g_second_half].value + beta_g_r1_first_half[g_first_half].value * beta_g_r1_second_half[g_second_half].value) % prime_field::mod * 
								(beta_u_first_half[u_first_half].value * beta_u_second_half[u_second_half].value % prime_field::mod) % prime_field::mod * (beta_v_first_half[v_first_half].value * beta_v_second_half[v_second_half].value % prime_field::mod);
					ret[1].value = ret[1].value % prime_field::mod;
					break;
				}
				case 2:
				{
					break;
				}
				case 3:
				{
					break;
				}
				case 4:
				{
					break;
				}
				case 5:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					
					auto beta_g_val = beta_g_r0_first_half[g_first_half] * beta_g_r0_second_half[g_second_half] + beta_g_r1_first_half[g_first_half] * beta_g_r1_second_half[g_second_half];
					auto beta_v_0 = beta_v_first_half[0] * beta_v_second_half[0];
					for(int j = u; j < v; ++j)
					{
						int u_first_half = j & ((1 << first_half_uv) - 1);
						int u_second_half = j >> first_half_uv;
						ret[5] = ret[5] + beta_g_val * beta_v_0 * (beta_u_first_half[u_first_half] * beta_u_second_half[u_second_half]);
					}
					break;
				}
				case 12:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					
					auto beta_g_val = beta_g_r0_first_half[g_first_half] * beta_g_r0_second_half[g_second_half] + beta_g_r1_first_half[g_first_half] * beta_g_r1_second_half[g_second_half];
					auto beta_v_0 = beta_v_first_half[0] * beta_v_second_half[0];
					for(int j = u; j <= v; ++j)
					{
						int u_first_half = j & ((1 << first_half_uv) - 1);
						int u_second_half = j >> first_half_uv;
						ret[12] = ret[12] + beta_g_val * beta_v_0 * (beta_u_first_half[u_first_half] * beta_u_second_half[u_second_half]);
						beta_v_0 = beta_v_0 + beta_v_0;
						if(beta_v_0.value >= prime_field::mod)
						{
							beta_v_0.value = beta_v_0.value - prime_field::mod;
						}
					}
					break;
				}
				case 6:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					int u_first_half = u & ((1 << first_half_uv) - 1);
					int u_second_half = u >> first_half_uv;
					int v_first_half = v & ((1 << first_half_uv) - 1);
					int v_second_half = v >> first_half_uv;
					ret[6].value = ret[6].value + (beta_g_r0_first_half[g_first_half].value * beta_g_r0_second_half[g_second_half].value + beta_g_r1_first_half[g_first_half].value * beta_g_r1_second_half[g_second_half].value) % prime_field::mod * 
								(beta_u_first_half[u_first_half].value * beta_u_second_half[u_second_half].value % prime_field::mod) % prime_field::mod * (beta_v_first_half[v_first_half].value * beta_v_second_half[v_second_half].value % prime_field::mod);
					ret[6].value = ret[6].value % prime_field::mod;
					break;
				}
				case 7:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					int u_first_half = u & ((1 << first_half_uv) - 1);
					int u_second_half = u >> first_half_uv;
					int v_first_half = v & ((1 << first_half_uv) - 1);
					int v_second_half = v >> first_half_uv;
					ret[7].value = ret[7].value + (beta_g_r0_first_half[g_first_half].value * beta_g_r0_second_half[g_second_half].value + beta_g_r1_first_half[g_first_half].value * beta_g_r1_second_half[g_second_half].value) % prime_field::mod * 
								(beta_u_first_half[u_first_half].value * beta_u_second_half[u_second_half].value % prime_field::mod) % prime_field::mod * (beta_v_first_half[v_first_half].value * beta_v_second_half[v_second_half].value % prime_field::mod);
					ret[7].value = ret[7].value % prime_field::mod;
					break;
				}
				case 8:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					int u_first_half = u & ((1 << first_half_uv) - 1);
					int u_second_half = u >> first_half_uv;
					int v_first_half = v & ((1 << first_half_uv) - 1);
					int v_second_half = v >> first_half_uv;
					ret[8].value = ret[8].value + (beta_g_r0_first_half[g_first_half].value * beta_g_r0_second_half[g_second_half].value + beta_g_r1_first_half[g_first_half].value * beta_g_r1_second_half[g_second_half].value) % prime_field::mod * 
								(beta_u_first_half[u_first_half].value * beta_u_second_half[u_second_half].value % prime_field::mod) % prime_field::mod * (beta_v_first_half[v_first_half].value * beta_v_second_half[v_second_half].value % prime_field::mod);
					ret[8].value = ret[8].value % prime_field::mod;
					break;
				}
				case 9:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					int u_first_half = u & ((1 << first_half_uv) - 1);
					int u_second_half = u >> first_half_uv;
					int v_first_half = v & ((1 << first_half_uv) - 1);
					int v_second_half = v >> first_half_uv;
					ret[9].value = ret[9].value + (beta_g_r0_first_half[g_first_half].value * beta_g_r0_second_half[g_second_half].value + beta_g_r1_first_half[g_first_half].value * beta_g_r1_second_half[g_second_half].value) % prime_field::mod * 
								(beta_u_first_half[u_first_half].value * beta_u_second_half[u_second_half].value % prime_field::mod) % prime_field::mod * (beta_v_first_half[v_first_half].value * beta_v_second_half[v_second_half].value % prime_field::mod);
					ret[9].value = ret[9].value % prime_field::mod;
					break;
				}
				case 10:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					ret[10].value = ret[10].value + (beta_g_r0_first_half[g_first_half].value * beta_g_r0_second_half[g_second_half].value + beta_g_r1_first_half[g_first_half].value * beta_g_r1_second_half[g_second_half].value) % prime_field::mod * tmp_u_val[u].value;
					ret[10].value = ret[10].value % prime_field::mod;
					break;
				}
				case 13:
				{
					int g_first_half = g & ((1 << first_half_g) - 1);
					int g_second_half = (g >> first_half_g);
					int u_first_half = u & ((1 << first_half_uv) - 1);
					int u_second_half = u >> first_half_uv;
					int v_first_half = v & ((1 << first_half_uv) - 1);
					int v_second_half = v >> first_half_uv;
					ret[13].value = ret[13].value + (beta_g_r0_first_half[g_first_half].value * beta_g_r0_second_half[g_second_half].value + beta_g_r1_first_half[g_first_half].value * beta_g_r1_second_half[g_second_half].value) % prime_field::mod * 
								(beta_u_first_half[u_first_half].value * beta_u_second_half[u_second_half].value % prime_field::mod) % prime_field::mod * (beta_v_first_half[v_first_half].value * beta_v_second_half[v_second_half].value % prime_field::mod);
					ret[13].value = ret[13].value % prime_field::mod;
					break;
				}
			}
		}
	}
	for(int i = 0; i < gate_type_count; ++i)
	{
		ret_para[i].value = ret_para[i].value % prime_field::mod;
		ret[i].value = ret[i].value % prime_field::mod;
		if(ret[i].value < 0)
			ret[i].value = ret[i].value + prime_field::mod;
		if(ret_para[i].value < 0)
			ret_para[i].value = ret_para[i].value + prime_field::mod;
		if(C.circuit[depth].is_parallel)
			assert(ret[i] == ret_para[i]);
	}
	return ret;
}

prime_field::field_element zk_verifier::direct_relay(int depth, prime_field::field_element *r_g, prime_field::field_element *r_u)
{
	if(depth != 1)
		return prime_field::field_element(0);
	else
	{
		prime_field::field_element ret = prime_field::field_element(1);
		for(int i = 0; i < C.circuit[depth].bit_length; ++i)
			ret = ret * (prime_field::field_element(1) - r_g[i] - r_u[i] + prime_field::field_element(2) * r_g[i] * r_u[i]);
		return ret;
	}
}

void zk_verifier::beta_init(int depth, prime_field::field_element alpha, prime_field::field_element beta,
	const prime_field::field_element* r_0, const prime_field::field_element* r_1, 
	const prime_field::field_element* r_u, const prime_field::field_element* r_v, 
	const prime_field::field_element* one_minus_r_0, const prime_field::field_element* one_minus_r_1, 
	const prime_field::field_element* one_minus_r_u, const prime_field::field_element* one_minus_r_v)
{
	bool debug_mode = false;
	if(!C.circuit[depth].is_parallel || debug_mode)
	{
		beta_g_r0_first_half[0] = alpha;
		beta_g_r1_first_half[0] = beta;
		beta_g_r0_second_half[0] = prime_field::field_element(1);
		beta_g_r1_second_half[0] = prime_field::field_element(1);
		int first_half_len = C.circuit[depth].bit_length / 2;
		int second_half_len = C.circuit[depth].bit_length - first_half_len;
		for(int i = 0; i < first_half_len; ++i)
		{
			for(int j = 0; j < (1 << i); ++j)
			{
				beta_g_r0_first_half[j | (1 << i)].value = beta_g_r0_first_half[j].value * r_0[i].value % prime_field::mod;
				beta_g_r1_first_half[j | (1 << i)].value = beta_g_r1_first_half[j].value * r_1[i].value % prime_field::mod;
			}
			for(int j = 0; j < (1 << i); ++j)
			{
				beta_g_r0_first_half[j].value = beta_g_r0_first_half[j].value * one_minus_r_0[i].value % prime_field::mod;
				beta_g_r1_first_half[j].value = beta_g_r1_first_half[j].value * one_minus_r_1[i].value % prime_field::mod;
			}
		}
		for(int i = 0; i < second_half_len; ++i)
		{
			for(int j = 0; j < (1 << i); ++j)
			{
				beta_g_r0_second_half[j | (1 << i)].value = beta_g_r0_second_half[j].value * r_0[i + first_half_len].value % prime_field::mod;
				beta_g_r1_second_half[j | (1 << i)].value = beta_g_r1_second_half[j].value * r_1[i + first_half_len].value % prime_field::mod;
			}
			for(int j = 0; j < (1 << i); ++j)
			{
				beta_g_r0_second_half[j].value = beta_g_r0_second_half[j].value * one_minus_r_0[i + first_half_len].value % prime_field::mod;
				beta_g_r1_second_half[j].value = beta_g_r1_second_half[j].value * one_minus_r_1[i + first_half_len].value % prime_field::mod;
			}
		}

		beta_u_first_half[0] = prime_field::field_element(1);
		beta_v_first_half[0] = prime_field::field_element(1);
		beta_u_second_half[0] = prime_field::field_element(1);
		beta_v_second_half[0] = prime_field::field_element(1);
		first_half_len = C.circuit[depth - 1].bit_length / 2;
		second_half_len = C.circuit[depth - 1].bit_length - first_half_len;

		for(int i = 0; i < first_half_len; ++i)
		{
			for(int j = 0; j < (1 << i); ++j)
			{
				beta_u_first_half[j | (1 << i)] = beta_u_first_half[j] * r_u[i];
				beta_v_first_half[j | (1 << i)] = beta_v_first_half[j] * r_v[i];
			}
				
			for(int j = 0; j < (1 << i); ++j)
			{
				beta_u_first_half[j] = beta_u_first_half[j] * one_minus_r_u[i];
				beta_v_first_half[j] = beta_v_first_half[j] * one_minus_r_v[i];
			}
		}

		for(int i = 0; i < second_half_len; ++i)
		{
			for(int j = 0; j < (1 << i); ++j)
			{
				beta_u_second_half[j | (1 << i)] = beta_u_second_half[j] * r_u[i + first_half_len];
				beta_v_second_half[j | (1 << i)] = beta_v_second_half[j] * r_v[i + first_half_len];
			}
				
			for(int j = 0; j < (1 << i); ++j)
			{
				beta_u_second_half[j] = beta_u_second_half[j] * one_minus_r_u[i + first_half_len];
				beta_v_second_half[j] = beta_v_second_half[j] * one_minus_r_v[i + first_half_len];
			}
		}
	}
	if(C.circuit[depth].is_parallel)
	{
		beta_g_r0_block_first_half[0] = alpha;
		beta_g_r1_block_first_half[0] = beta;
		beta_g_r0_block_second_half[0] = prime_field::field_element(1);
		beta_g_r1_block_second_half[0] = prime_field::field_element(1);
		int first_half_len = C.circuit[depth].log_block_size / 2;
		int second_half_len = C.circuit[depth].log_block_size - first_half_len;
		for(int i = 0; i < first_half_len; ++i)
		{
			for(int j = 0; j < (1 << i); ++j)
			{
				beta_g_r0_block_first_half[j | (1 << i)].value = beta_g_r0_block_first_half[j].value * r_0[i].value % prime_field::mod;
				beta_g_r1_block_first_half[j | (1 << i)].value = beta_g_r1_block_first_half[j].value * r_1[i].value % prime_field::mod;
			}
			for(int j = 0; j < (1 << i); ++j)
			{
				beta_g_r0_block_first_half[j].value = beta_g_r0_block_first_half[j].value * one_minus_r_0[i].value % prime_field::mod;
				beta_g_r1_block_first_half[j].value = beta_g_r1_block_first_half[j].value * one_minus_r_1[i].value % prime_field::mod;
			}
		}
		for(int i = 0; i < second_half_len; ++i)
		{
			for(int j = 0; j < (1 << i); ++j)
			{
				beta_g_r0_block_second_half[j | (1 << i)].value = beta_g_r0_block_second_half[j].value * r_0[i + first_half_len].value % prime_field::mod;
				beta_g_r1_block_second_half[j | (1 << i)].value = beta_g_r1_block_second_half[j].value * r_1[i + first_half_len].value % prime_field::mod;
			}
			for(int j = 0; j < (1 << i); ++j)
			{
				beta_g_r0_block_second_half[j].value = beta_g_r0_block_second_half[j].value * one_minus_r_0[i + first_half_len].value % prime_field::mod;
				beta_g_r1_block_second_half[j].value = beta_g_r1_block_second_half[j].value * one_minus_r_1[i + first_half_len].value % prime_field::mod;
			}
		}

		beta_u_block_first_half[0] = prime_field::field_element(1);
		beta_v_block_first_half[0] = prime_field::field_element(1);
		beta_u_block_second_half[0] = prime_field::field_element(1);
		beta_v_block_second_half[0] = prime_field::field_element(1);
		first_half_len = C.circuit[depth - 1].log_block_size / 2;
		second_half_len = C.circuit[depth - 1].log_block_size - first_half_len;

		for(int i = 0; i < first_half_len; ++i)
		{
			for(int j = 0; j < (1 << i); ++j)
			{
				beta_u_block_first_half[j | (1 << i)] = beta_u_block_first_half[j] * r_u[i];
				beta_v_block_first_half[j | (1 << i)] = beta_v_block_first_half[j] * r_v[i];
			}
				
			for(int j = 0; j < (1 << i); ++j)
			{
				beta_u_block_first_half[j] = beta_u_block_first_half[j] * one_minus_r_u[i];
				beta_v_block_first_half[j] = beta_v_block_first_half[j] * one_minus_r_v[i];
			}
		}

		for(int i = 0; i < second_half_len; ++i)
		{
			for(int j = 0; j < (1 << i); ++j)
			{
				beta_u_block_second_half[j | (1 << i)] = beta_u_block_second_half[j] * r_u[i + first_half_len];
				beta_v_block_second_half[j | (1 << i)] = beta_v_block_second_half[j] * r_v[i + first_half_len];
			}
				
			for(int j = 0; j < (1 << i); ++j)
			{
				beta_u_block_second_half[j] = beta_u_block_second_half[j] * one_minus_r_u[i + first_half_len];
				beta_v_block_second_half[j] = beta_v_block_second_half[j] * one_minus_r_v[i + first_half_len];
			}
		}
	}
}

prime_field::field_element* zk_verifier::generate_randomness(unsigned int size)
{
	int k = size;
	prime_field::field_element* ret;
	ret = trans.random(size);
	return ret;
}

prime_field::field_element zk_verifier::V_in(const prime_field::field_element* r_0, const prime_field::field_element* one_minus_r_0,
								prime_field::field_element* output_raw, int r_0_size, int output_size)
{
	prime_field::field_element* output = new prime_field::field_element[output_size];
	for(int i = 0; i < output_size; ++i)
		output[i] = output_raw[i];
	for(int i = 0; i < r_0_size; ++i)
	{
		for(int j = 0; j < (output_size >> 1); ++j)
			output[j] = output[j << 1] * (one_minus_r_0[i]) + output[j << 1 | 1] * (r_0[i]);
		output_size >>= 1;
	}
	auto ret = output[0];
	ret.value = ret.value % prime_field::mod;
	delete[] output;
	if(ret.value < 0)
		ret.value = ret.value + prime_field::mod;
	return ret;
}

bool zk_verifier::verify(const char* output_path)
{
	int proof_size = 0;
	//there is a way to compress binlinear pairing element
	int bilinear_pairing_factor = 3;
	double verification_time = 0;
	double predicates_calc_time = 0;
	double verification_rdl_time = 0;
	double poly_commit_scheme_time = 0;
	prime_field::init_random();
	p -> proof_init();

	auto result = p -> evaluate();
	double key_gen_time = 0;
	auto commit_t0 = std::chrono::high_resolution_clock::now();
	auto digest_input = p -> keygen_and_commit(C.circuit[0].bit_length, key_gen_time);
	auto commit_t1 = std::chrono::high_resolution_clock::now();
	poly_commit_scheme_time += (std::chrono::duration_cast<std::chrono::duration<double>>(commit_t1 - commit_t0)).count() - key_gen_time;
	proof_size += sizeof(bn::Ec1) / bilinear_pairing_factor * (digest_input.first.size() + digest_input.second.size());
	for(int i = 0; i < digest_input.first.size(); ++i)
		trans.msg.push_back(container(digest_input.first[i]));
	for(int i = 0; i < digest_input.second.size(); ++i)
		trans.msg.push_back(container(digest_input.second[i]));
	prime_field::field_element alpha, beta;
	alpha.value = 1;
	beta.value = 0;
	random_oracle oracle;
	//initial random value
	prime_field::field_element *r_0 = generate_randomness(C.circuit[C.total_depth - 1].bit_length), *r_1 = generate_randomness(C.circuit[C.total_depth - 1].bit_length);
	prime_field::field_element *one_minus_r_0, *one_minus_r_1;
	one_minus_r_0 = new prime_field::field_element[C.circuit[C.total_depth - 1].bit_length];
	one_minus_r_1 = new prime_field::field_element[C.circuit[C.total_depth - 1].bit_length];

	for(int i = 0; i < (C.circuit[C.total_depth - 1].bit_length); ++i)
	{
		one_minus_r_0[i] = prime_field::field_element(1) - r_0[i];
		one_minus_r_1[i] = prime_field::field_element(1) - r_1[i];
	}
	
	std::chrono::high_resolution_clock::time_point t_a = std::chrono::high_resolution_clock::now();
	std::cerr << "Calc V_output(r)" << std::endl;
	prime_field::field_element a_0 = p -> V_res(one_minus_r_0, r_0, result, C.circuit[C.total_depth - 1].bit_length, (1 << (C.circuit[C.total_depth - 1].bit_length)));
	std::chrono::high_resolution_clock::time_point t_b = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> ts = std::chrono::duration_cast<std::chrono::duration<double>>(t_b - t_a);
	std::cerr << "	Time: " << ts.count() << std::endl;
	a_0 = alpha * a_0;

	prime_field::field_element alpha_beta_sum = a_0; //+ a_1

	prime_field::field_element direct_relay_value;
	for(int i = C.total_depth - 1; i >= 1; --i)
	{
		std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
	//	std::cerr << "Bound u start" << std::endl;
		auto rho = trans.random();
		std::vector<bn::Ec1> digest_mask;

		auto digest_maskR = p -> sumcheck_init(i, C.circuit[i].bit_length, C.circuit[i - 1].bit_length, C.circuit[i - 1].bit_length, alpha, beta, r_0, r_1, one_minus_r_0, one_minus_r_1);

		digest_mask = p -> generate_maskpoly_pre_rho(C.circuit[i - 1].bit_length * 2 + 1, 2);
		proof_size += sizeof(bn::Ec1) / bilinear_pairing_factor * (digest_mask.size() + digest_maskR.size());
		for(int j = 0; j < digest_mask.size(); ++j)
			trans.msg.push_back(digest_mask[j]);
		for(int j = 0; j < digest_maskR.size(); ++j)
			trans.msg.push_back(digest_maskR[j]);
		p -> rho = rho;
		p -> generate_maskpoly_after_rho(C.circuit[i - 1].bit_length * 2 + 1, 2);
		bool r_verify_cc = vpdR::check_commit(digest_maskR[0], digest_maskR[1]);
		bool msk_poly_cc = vpd_test::check_commit(digest_mask[0], digest_mask[1]);
		
		//add maskpoly

		alpha_beta_sum.value = (alpha_beta_sum.value + p->maskpoly_sumc.value) % prime_field::mod;

		p -> sumcheck_phase1_init();
		prime_field::field_element previous_random = prime_field::field_element(0);
		//next level random
		auto r_u = generate_randomness(C.circuit[i - 1].bit_length);
		auto r_v = generate_randomness(C.circuit[i - 1].bit_length);
		direct_relay_value = alpha * direct_relay(i, r_0, r_u) + beta * direct_relay(i, r_1, r_u);
		auto r_c = generate_randomness(1);
		prime_field::field_element *one_minus_r_u, *one_minus_r_v;
		one_minus_r_u = new prime_field::field_element[C.circuit[i - 1].bit_length];
		one_minus_r_v = new prime_field::field_element[C.circuit[i - 1].bit_length];
		
		for(int j = 0; j < C.circuit[i - 1].bit_length; ++j)
		{
			one_minus_r_u[j] = prime_field::field_element(1) - r_u[j];
			one_minus_r_v[j] = prime_field::field_element(1) - r_v[j];
		}

		for(int j = 0; j < C.circuit[i - 1].bit_length; ++j)
		{	
			if(j == C.circuit[i - 1].bit_length - 1){
				quintuple_poly poly = p->sumcheck_phase1_updatelastbit(previous_random, j);
				proof_size += sizeof(quintuple_poly) / 2;
				trans.msg.push_back(container(poly.a, container::fld_ele_indicator));
				trans.msg.push_back(container(poly.b, container::fld_ele_indicator));
				trans.msg.push_back(container(poly.c, container::fld_ele_indicator));
				trans.msg.push_back(container(poly.d, container::fld_ele_indicator));
				trans.msg.push_back(container(poly.e, container::fld_ele_indicator));
				trans.msg.push_back(container(poly.f, container::fld_ele_indicator));
				previous_random = r_u[j];


				if(poly.eval(0) + poly.eval(1) != alpha_beta_sum)
				{ 
					std::cout << "round j = " << j << std::endl;
					fprintf(stderr, "Verification fail, phase1, circuit %d, current bit %d\n", i, j);
					return false;
				}
				else
				{
					
				}
				alpha_beta_sum = poly.eval(r_u[j]);
			}

			else{
				quadratic_poly poly = p -> sumcheck_phase1_update(previous_random, j);
				proof_size += sizeof(quadratic_poly) / 2;
				trans.msg.push_back(container(poly.a, container::fld_ele_indicator));
				trans.msg.push_back(container(poly.b, container::fld_ele_indicator));
				trans.msg.push_back(container(poly.c, container::fld_ele_indicator));
				previous_random = r_u[j];
			

				if(poly.eval(0) + poly.eval(1) != alpha_beta_sum)
				{
					std::cout << "round j = " << j << std::endl;
					fprintf(stderr, "Verification fail, phase1, circuit %d, current bit %d\n", i, j);
					return false;
				}
				else
				{
					
				}
				alpha_beta_sum = poly.eval(r_u[j]);
			}
		}
		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	//	std::cerr << "Bound v start" << std::endl;
		t0 = std::chrono::high_resolution_clock::now();
		p -> sumcheck_phase2_init(previous_random, r_u, one_minus_r_u);
		previous_random = prime_field::field_element(0);
		for(int j = 0; j < C.circuit[i - 1].bit_length; ++j)
		{
			if(i == 1)
				r_v[j] = prime_field::field_element(0);
			if(j == C.circuit[i - 1].bit_length - 1){
				quintuple_poly poly = p -> sumcheck_phase2_updatelastbit(previous_random, j);
				proof_size += sizeof(quintuple_poly) / 2;
				trans.msg.push_back(container(poly.a, container::fld_ele_indicator));
				trans.msg.push_back(container(poly.b, container::fld_ele_indicator));
				trans.msg.push_back(container(poly.c, container::fld_ele_indicator));
				trans.msg.push_back(container(poly.d, container::fld_ele_indicator));
				trans.msg.push_back(container(poly.e, container::fld_ele_indicator));
				trans.msg.push_back(container(poly.f, container::fld_ele_indicator));
				poly.f = poly.f;
				previous_random = r_v[j];
				if(poly.eval(0) + poly.eval(1) + direct_relay_value * p -> v_u != alpha_beta_sum)
				{
					fprintf(stderr, "Verification fail, phase2, circuit level %d, current bit %d\n", i, j);
					return false;
				}
				else
				{
					
				}
				alpha_beta_sum = poly.eval(r_v[j]) + direct_relay_value * p -> v_u;
			}
			else
			{
				quadratic_poly poly = p -> sumcheck_phase2_update(previous_random, j);
				proof_size += sizeof(quadratic_poly) / 2;
				trans.msg.push_back(container(poly.a, container::fld_ele_indicator));
				trans.msg.push_back(container(poly.b, container::fld_ele_indicator));
				trans.msg.push_back(container(poly.c, container::fld_ele_indicator));
				poly.c = poly.c;
			
				previous_random = r_v[j];
				if(poly.eval(0) + poly.eval(1) + direct_relay_value * p -> v_u != alpha_beta_sum)
				{
					fprintf(stderr, "Verification fail, phase2, circuit level %d, current bit %d\n", i, j);
					return false;
				}
				else
				{
					
				}
				alpha_beta_sum = poly.eval(r_v[j]) + direct_relay_value * p -> v_u;
			}
		}

	
		//Add one more round for maskR
		//quadratic_poly poly p->sumcheck_finalroundR(previous_random, C.current[i - 1].bit_length);
		
		auto final_claims = p -> sumcheck_finalize(previous_random);
		

		auto v_u = final_claims.first;
		auto v_v = final_claims.second;

		std::chrono::high_resolution_clock::time_point predicates_calc_t0 = std::chrono::high_resolution_clock::now();
		beta_init(i, alpha, beta, r_0, r_1, r_u, r_v, one_minus_r_0, one_minus_r_1, one_minus_r_u, one_minus_r_v);
		auto predicates_value = predicates(i, r_0, r_1, r_u, r_v, alpha, beta);
		std::chrono::high_resolution_clock::time_point predicates_calc_t1 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> predicates_calc_span = std::chrono::duration_cast<std::chrono::duration<double>>(predicates_calc_t1 - predicates_calc_t0);
		if(C.circuit[i].is_parallel == false)
			verification_rdl_time += predicates_calc_span.count();
		verification_time += predicates_calc_span.count();
		predicates_calc_time += predicates_calc_span.count();
		
		auto mult_value = predicates_value[1];
		auto add_value = predicates_value[0];
		auto not_value = predicates_value[6];
		auto minus_value = predicates_value[7];
		auto xor_value = predicates_value[8];
		auto naab_value = predicates_value[9];
		auto sum_value = predicates_value[5];
		auto relay_value = predicates_value[10];
		auto exp_sum_value = predicates_value[12];
		auto bit_test_value = predicates_value[13];
		quadratic_poly poly = p->sumcheck_finalround(previous_random, C.circuit[i - 1].bit_length << 1, add_value * (v_u + v_v) + mult_value * v_u * v_v + not_value * (prime_field::field_element(1) - v_u) + minus_value * (v_u - v_v) + xor_value * (v_u + v_v - prime_field::field_element(2) * v_u * v_v) + naab_value * (v_v - v_u * v_v) + sum_value * v_u + relay_value * v_u + exp_sum_value * v_u + bit_test_value * (v_u * (prime_field::field_element(1) - v_v)));

		if(poly.eval(0) + poly.eval(1) + direct_relay_value * v_u != alpha_beta_sum)
		{
			fprintf(stderr, "Verification fail, phase2, lastbit for c\n");
			return false;
		}
		if(i == 1)
			r_c[0] = prime_field::field_element(0);
		alpha_beta_sum = poly.eval(r_c[0]) + direct_relay_value * p -> v_u;

		auto poly_commit_prove_time_t0 = std::chrono::high_resolution_clock::now();
		mpz_class maskRg1_value_mpz, maskRg2_value_mpz;
		std::vector<mpz_class> r;
		r.resize(2);
		r[0] = p -> prepreu1.to_gmp_class(), r[1] = r_c[0].to_gmp_class();
		auto witnesses = p -> prove_R(r, maskRg1_value_mpz);
		prime_field::field_element tmp_rg1;
		proof_size += sizeof(bn::Ec1) / bilinear_pairing_factor * (witnesses.first.size() + witnesses.second.size());
		for(int i = 0; i < witnesses.first.size(); ++i)
			trans.msg.push_back(container(witnesses.first[i]));
		for(int i = 0; i < witnesses.second.size(); ++i)
			trans.msg.push_back(container(witnesses.second[i]));

		std::chrono::high_resolution_clock::time_point vpdr_verify_0_0 = std::chrono::high_resolution_clock::now();
		bool r_verify_verify = vpdR::verify(r, digest_maskR[0], maskRg1_value_mpz, witnesses.first, witnesses.second);
		std::chrono::high_resolution_clock::time_point vpdr_verify_0_1 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> vpdr_verify_0_span = std::chrono::duration_cast<std::chrono::duration<double>>(vpdr_verify_0_1 - vpdr_verify_0_0);
		verification_time += vpdr_verify_0_span.count();

		r[0] = p -> preprev1.to_gmp_class();
		witnesses = p -> prove_R(r, maskRg2_value_mpz);
		proof_size += sizeof(bn::Ec1) / bilinear_pairing_factor * (witnesses.first.size() + witnesses.second.size());
		for(int i = 0; i < witnesses.first.size(); ++i)
			trans.msg.push_back(container(witnesses.first[i]));
		for(int i = 0; i < witnesses.second.size(); ++i)
			trans.msg.push_back(container(witnesses.second[i]));
		
		std::chrono::high_resolution_clock::time_point vpdr_verify_1_0 = std::chrono::high_resolution_clock::now();
		r_verify_verify &= vpdR::verify(r, digest_maskR[0], maskRg2_value_mpz, witnesses.first, witnesses.second);
		std::chrono::high_resolution_clock::time_point vpdr_verify_1_1 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> vpdr_verify_1_span = std::chrono::duration_cast<std::chrono::duration<double>>(vpdr_verify_1_1 - vpdr_verify_1_0);
		verification_time += vpdr_verify_1_span.count();
		auto poly_commit_prove_time_t1 = std::chrono::high_resolution_clock::now();
		poly_commit_scheme_time += std::chrono::duration_cast<std::chrono::duration<double>>(poly_commit_prove_time_t1 - poly_commit_prove_time_t0).count();
		if(r_verify_verify & r_verify_cc)
		{
		}
		else
		{
			fprintf(stderr, "VPD R failed\n");
			return false;
		}

		r.clear();
		for(int j = 0; j < C.circuit[i - 1].bit_length; ++j)
			r.push_back(r_u[j].to_gmp_class());
		for(int j = 0; j < C.circuit[i - 1].bit_length; ++j)
			r.push_back(r_v[j].to_gmp_class());
		r.push_back(r_c[0].to_gmp_class());
		mpz_class maskpoly_value_mpz = 0;
		witnesses = p -> prove_mask(r, maskpoly_value_mpz);
		proof_size += sizeof(bn::Ec1) / bilinear_pairing_factor * (witnesses.first.size() + witnesses.second.size());
		for(int i = 0; i < witnesses.first.size(); ++i)
			trans.msg.push_back(container(witnesses.first[i]));
		for(int i = 0; i < witnesses.second.size(); ++i)
			trans.msg.push_back(container(witnesses.second[i]));



		std::chrono::high_resolution_clock::time_point vpd_test_verify_0 = std::chrono::high_resolution_clock::now();
		auto msk_poly_verify = vpd_test::verify(r, digest_mask[0], maskpoly_value_mpz, witnesses.first, witnesses.second);
		std::chrono::high_resolution_clock::time_point vpd_test_verify_1 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> vpd_test_verify_span = std::chrono::duration_cast<std::chrono::duration<double>>(vpd_test_verify_1 - vpd_test_verify_0);
		verification_time += vpd_test_verify_span.count();

		if(msk_poly_cc & msk_poly_verify)
		{
		}
		else
		{
			fprintf(stderr, "VPD mask fail\n");
			return false;
		}
		
		prime_field::field_element maskpoly_value;
		prime_field::field_element maskRg1_value;
		prime_field::field_element maskRg2_value;
		maskpoly_value.value = prime_field::u512b(maskpoly_value_mpz.get_str().c_str(), maskpoly_value_mpz.get_str().length(), 10);
		maskRg1_value.value = prime_field::u512b(maskRg1_value_mpz.get_str().c_str(), maskRg1_value_mpz.get_str().length(), 10);
		maskRg2_value.value = prime_field::u512b(maskRg2_value_mpz.get_str().c_str(), maskRg2_value_mpz.get_str().length(), 10);
		if(alpha_beta_sum != r_c[0] * (add_value * (v_u + v_v) + mult_value * v_u * v_v + not_value * (prime_field::field_element(1) - v_u) + minus_value * (v_u - v_v) + xor_value * (v_u + v_v - prime_field::field_element(2) * v_u * v_v) + naab_value * (v_v - v_u * v_v) + sum_value * v_u + relay_value * v_u + exp_sum_value * v_u + bit_test_value * (prime_field::field_element(1) - v_v) * v_u) + alpha * p -> Iuv * p ->preZu * maskRg1_value + beta * p -> Iuv * p -> preZv * maskRg2_value + rho * maskpoly_value + direct_relay_value * v_u)
		{
			fprintf(stderr, "Verification fail, semi final, circuit level %d\n", i);
			return false;
		}
		else
		{
		}
		auto tmp_alpha = generate_randomness(1), tmp_beta = generate_randomness(1);
		alpha = tmp_alpha[0];
		beta = tmp_beta[0];
		delete[] tmp_alpha;
		delete[] tmp_beta;
		if(i != 1)
			alpha_beta_sum = alpha * v_u + beta * v_v;
		else
			alpha_beta_sum = v_u;
		delete[] r_0;
		delete[] r_1;
		delete[] one_minus_r_0;
		delete[] one_minus_r_1;
		r_0 = r_u;
		r_1 = r_v;
		one_minus_r_0 = one_minus_r_u;
		one_minus_r_1 = one_minus_r_v;
	}

	//post sumcheck

	prime_field::field_element input_0;
	
	std::vector<mpz_class> r_0_mpz, r_1_mpz;
	for(int i = 0; i< C.circuit[0].bit_length; ++i)
		r_0_mpz.push_back(r_0[i].to_gmp_class());
	for(int i = 0; i< C.circuit[0].bit_length; ++i)
		r_1_mpz.push_back(r_1[i].to_gmp_class());
	
	mpz_class input_0_mpz, input_1_mpz;
	auto poly_commit_prove_time_t0 = std::chrono::high_resolution_clock::now();

	input_0_mpz = 0, input_1_mpz = 0;
	auto witnesses_0 = p -> prove_input(r_0_mpz, input_0_mpz, p -> Zu.to_gmp_class());
	proof_size += sizeof(bn::Ec1) / bilinear_pairing_factor * (witnesses_0.first.size() + witnesses_0.second.size());
	for(int i = 0; i < witnesses_0.first.size(); ++i)
		trans.msg.push_back(container(witnesses_0.first[i]));
	for(int i = 0; i < witnesses_0.second.size(); ++i)
		trans.msg.push_back(container(witnesses_0.second[i]));
	

	std::chrono::high_resolution_clock::time_point vpd_input_0 = std::chrono::high_resolution_clock::now();
	bool input_0_verify = input_vpd::verify(r_0_mpz, digest_input.first[0], digest_input.second[0], p -> Zu.to_gmp_class(), input_0_mpz, witnesses_0.first, witnesses_0.second);
	std::chrono::high_resolution_clock::time_point vpd_input_1 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> vpd_input_span = std::chrono::duration_cast<std::chrono::duration<double>>(vpd_input_1 - vpd_input_0);
	verification_time += vpd_input_span.count();
	auto poly_commit_prove_time_t1 = std::chrono::high_resolution_clock::now();
	poly_commit_scheme_time += std::chrono::duration_cast<std::chrono::duration<double>>(poly_commit_prove_time_t1 - poly_commit_prove_time_t0).count();
	if(!(input_0_verify))
	{
		fprintf(stderr, "Verification fail, input vpd.\n");
		return false;
	}

	input_0 = input_0 + p->Zu * p->sumRc.eval(p->preu1);
	
	auto is0 = input_0_mpz.get_str();
	input_0.value = prime_field::u512b(is0.c_str(), is0.length(), 10);
	
	delete[] r_0;
	delete[] r_1;
	delete[] one_minus_r_0;
	delete[] one_minus_r_1;
	if(alpha_beta_sum != input_0)
	{
		fprintf(stderr, "Verification fail, final input check fail.\n");
		return false;
	}
	else
	{
		fprintf(stderr, "Verification pass\n");
		std::cerr << "Prove Time " << p -> total_time << std::endl;
	//	std::cerr << "Verification gate time " << predicates_calc_time << std::endl;
		std::cerr << "Verification rdl time " << verification_rdl_time << std::endl;
		std::cerr << "Verification Time " << verification_time - verification_rdl_time << std::endl;
		std::cerr << "Proof size(bytes) est " << proof_size << std::endl;
		std::cerr << "Proof size(bytes) real " << trans.get_proof_size() << std::endl;
		std::cerr << "Poly commit scheme time " << poly_commit_scheme_time << std::endl;
		trans.output_to_file("proof.bin");
		FILE *result = fopen(output_path, "w");
		fprintf(result, "%lf %lf %lf %lf %lf %d %lf\n", p -> total_time, verification_time, predicates_calc_time, verification_rdl_time, key_gen_time, proof_size, poly_commit_scheme_time);
		fclose(result);
	}
	p -> delete_self();
	delete_self();
	return true;
}

void zk_verifier::delete_self()
{
	delete[] beta_g_r0_first_half;
	delete[] beta_g_r0_second_half;
	delete[] beta_g_r1_first_half;
	delete[] beta_g_r1_second_half;
	delete[] beta_u_first_half;
	delete[] beta_u_second_half;
	delete[] beta_v_first_half;
	delete[] beta_v_second_half;
	delete[] beta_g_r0_block_first_half;
	delete[] beta_g_r0_block_second_half;
	delete[] beta_g_r1_block_first_half;
	delete[] beta_g_r1_block_second_half;
	delete[] beta_u_block_first_half;
	delete[] beta_u_block_second_half;
	delete[] beta_v_block_first_half;
	delete[] beta_v_block_second_half;
	for(int i = 0; i < C.total_depth; ++i)
	{
		delete[] C.circuit[i].gates;
	}
	delete[] C.circuit;
}
