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

void zk_verifier::read_circuit(const char *path)
{
	int d;
	FILE *circuit_in;
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
			if(ty == 6)
			{
				if(v != 0)
					fprintf(stderr, "WARNING, v!=0 for NOT gate.\n");
				v = 0;
			}
			if(g != previous_g + 1)
			{
				printf("Error, gates must be in sorted order, and full [0, 2^n - 1].");
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
	fclose(circuit_in);
}

vector<prime_field::field_element> zk_verifier::predicates(int depth)
{
	vector<prime_field::field_element> ret;
	ret.resize(10);
	for(int i = 0; i < 10; ++i)
		ret[i] = prime_field::field_element(0);
	int first_half_g = C.circuit[depth].bit_length / 2;
	int second_half_g = C.circuit[depth].bit_length - first_half_g;
	int first_half_uv = C.circuit[depth - 1].bit_length / 2;
	int second_half_uv = C.circuit[depth - 1].bit_length - first_half_uv;
	if(depth == 1)
	{
		return ret;
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
				ret[0] = ret[0] + (beta_g_r0_first_half[g_first_half] * beta_g_r0_second_half[g_second_half] + beta_g_r1_first_half[g_first_half] * beta_g_r1_second_half[g_second_half]) * 
							(beta_u_first_half[u_first_half] * beta_u_second_half[u_second_half]) * (beta_v_first_half[v_first_half] * beta_v_second_half[v_second_half]);
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
				ret[1] = ret[1] + (beta_g_r0_first_half[g_first_half] * beta_g_r0_second_half[g_second_half] + beta_g_r1_first_half[g_first_half] * beta_g_r1_second_half[g_second_half]) * 
							(beta_u_first_half[u_first_half] * beta_u_second_half[u_second_half]) * (beta_v_first_half[v_first_half] * beta_v_second_half[v_second_half]);
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
			case 6:
			{
				int g_first_half = g & ((1 << first_half_g) - 1);
				int g_second_half = (g >> first_half_g);
				int u_first_half = u & ((1 << first_half_uv) - 1);
				int u_second_half = u >> first_half_uv;
				int v_first_half = v & ((1 << first_half_uv) - 1);
				int v_second_half = v >> first_half_uv;
				ret[6] = ret[6] + (beta_g_r0_first_half[g_first_half] * beta_g_r0_second_half[g_second_half] + beta_g_r1_first_half[g_first_half] * beta_g_r1_second_half[g_second_half]) * 
							(beta_u_first_half[u_first_half] * beta_u_second_half[u_second_half]) * (beta_v_first_half[v_first_half] * beta_v_second_half[v_second_half]);
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
				ret[7] = ret[7] + (beta_g_r0_first_half[g_first_half] * beta_g_r0_second_half[g_second_half] + beta_g_r1_first_half[g_first_half] * beta_g_r1_second_half[g_second_half]) * 
							(beta_u_first_half[u_first_half] * beta_u_second_half[u_second_half]) * (beta_v_first_half[v_first_half] * beta_v_second_half[v_second_half]);
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
				ret[8] = ret[8] + (beta_g_r0_first_half[g_first_half] * beta_g_r0_second_half[g_second_half] + beta_g_r1_first_half[g_first_half] * beta_g_r1_second_half[g_second_half]) * 
							(beta_u_first_half[u_first_half] * beta_u_second_half[u_second_half]) * (beta_v_first_half[v_first_half] * beta_v_second_half[v_second_half]);
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
				ret[9] = ret[9] + (beta_g_r0_first_half[g_first_half] * beta_g_r0_second_half[g_second_half] + beta_g_r1_first_half[g_first_half] * beta_g_r1_second_half[g_second_half]) * 
							(beta_u_first_half[u_first_half] * beta_u_second_half[u_second_half]) * (beta_v_first_half[v_first_half] * beta_v_second_half[v_second_half]);
				break;
			}
		}
	}
	for(int i = 0; i < 10; ++i)
	{
		ret[i].value = ret[i].value % prime_field::mod;
		if(ret[i].value < 0)
			ret[i].value = ret[i].value + prime_field::mod;
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

prime_field::field_element* generate_randomness(unsigned int size)
{
	int k = size;
	prime_field::field_element* ret;
	ret = new prime_field::field_element[k];

	for(int i = 0; i < k; ++i)
	{
		ret[i] = prime_field::random();
		ret[i].value = ret[i].value % prime_field::mod;
	}
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

bool zk_verifier::verify()
{
	double verification_time = 0;
	prime_field::init_random();
	p -> proof_init();

	auto result = p -> evaluate();

	auto digest_input = p -> keygen_and_commit(C.circuit[0].bit_length);

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
		std::cerr << "Bound u start" << std::endl;
		auto rho = prime_field::random();
		std::vector<bn::Ec1> digest_mask;

		auto digest_maskR = p -> sumcheck_init(i, C.circuit[i].bit_length, C.circuit[i - 1].bit_length, C.circuit[i - 1].bit_length, alpha, beta, r_0, r_1, one_minus_r_0, one_minus_r_1);
		digest_mask = p -> generate_maskpoly_pre_rho(C.circuit[i - 1].bit_length * 2 + 1, 2);
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
		std::cerr << "	Time: " << time_span.count() << std::endl;
		std::cerr << "Bound v start" << std::endl;
		t0 = std::chrono::high_resolution_clock::now();
		p -> sumcheck_phase2_init(previous_random, r_u, one_minus_r_u);
		previous_random = prime_field::field_element(0);
		for(int j = 0; j < C.circuit[i - 1].bit_length; ++j)
		{
			if(i == 1)
				r_v[j] = prime_field::field_element(0);
			if(j == C.circuit[i - 1].bit_length - 1){
				quintuple_poly poly = p -> sumcheck_phase2_updatelastbit(previous_random, j);
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
		auto predicates_value = predicates(i);
		std::chrono::high_resolution_clock::time_point predicates_calc_t1 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> predicates_calc_span = std::chrono::duration_cast<std::chrono::duration<double>>(predicates_calc_t1 - predicates_calc_t0);
		verification_time += predicates_calc_span.count();
		
		auto mult_value = predicates_value[1];
		auto add_value = predicates_value[0];
		auto not_value = predicates_value[6];
		auto minus_value = predicates_value[7];
		auto xor_value = predicates_value[8];
		auto naab_value = predicates_value[9];
		auto sum_value = predicates_value[5];
		quadratic_poly poly = p->sumcheck_finalround(previous_random, C.circuit[i - 1].bit_length << 1, add_value * (v_u + v_v) + mult_value * v_u * v_v + not_value * (prime_field::field_element(1) - v_u) + minus_value * (v_u - v_v) + xor_value * (v_u + v_v - prime_field::field_element(2) * v_u * v_v) + naab_value * (v_v - v_u * v_v) + sum_value * v_u);

		if(poly.eval(0) + poly.eval(1) + direct_relay_value * v_u != alpha_beta_sum)
		{
			fprintf(stderr, "Verification fail, phase2, lastbit for c\n");
			return false;
		}
		if(i == 1)
			r_c[0] = prime_field::field_element(0);
		alpha_beta_sum = poly.eval(r_c[0]) + direct_relay_value * p -> v_u;

		mpz_class maskRg1_value_mpz, maskRg2_value_mpz;
		std::vector<mpz_class> r;
		r.resize(2);
		r[0] = p -> prepreu1.to_gmp_class(), r[1] = r_c[0].to_gmp_class();
		auto witnesses = p -> prove_R(r, maskRg1_value_mpz);
		prime_field::field_element tmp_rg1;


		std::chrono::high_resolution_clock::time_point vpdr_verify_0_0 = std::chrono::high_resolution_clock::now();
		bool r_verify_verify = vpdR::verify(r, digest_maskR[0], maskRg1_value_mpz, witnesses.first, witnesses.second);
		std::chrono::high_resolution_clock::time_point vpdr_verify_0_1 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> vpdr_verify_0_span = std::chrono::duration_cast<std::chrono::duration<double>>(vpdr_verify_0_1 - vpdr_verify_0_0);
		verification_time += vpdr_verify_0_span.count();

		r[0] = p -> preprev1.to_gmp_class();
		witnesses = p -> prove_R(r, maskRg2_value_mpz);
		
		std::chrono::high_resolution_clock::time_point vpdr_verify_1_0 = std::chrono::high_resolution_clock::now();
		r_verify_verify &= vpdR::verify(r, digest_maskR[0], maskRg2_value_mpz, witnesses.first, witnesses.second);
		std::chrono::high_resolution_clock::time_point vpdr_verify_1_1 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> vpdr_verify_1_span = std::chrono::duration_cast<std::chrono::duration<double>>(vpdr_verify_1_1 - vpdr_verify_1_0);
		verification_time += vpdr_verify_1_span.count();

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
		if(alpha_beta_sum != r_c[0] * (add_value * (v_u + v_v) + mult_value * v_u * v_v + not_value * (prime_field::field_element(1) - v_u) + minus_value * (v_u - v_v) + xor_value * (v_u + v_v - prime_field::field_element(2) * v_u * v_v) + naab_value * (v_v - v_u * v_v) + sum_value * v_u) + alpha * p -> Iuv * p ->preZu * maskRg1_value + beta * p -> Iuv * p -> preZv * maskRg2_value + rho * maskpoly_value + direct_relay_value * v_u)
		{
			fprintf(stderr, "Verification fail, semi final, circuit level %d\n", i);
			return false;
		}
		else
		{
			fprintf(stderr, "Verification Pass, semi final, circuit level %d\n", i);
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

	input_0_mpz = 0, input_1_mpz = 0;
	auto witnesses_0 = p -> prove_input(r_0_mpz, input_0_mpz, p -> Zu.to_gmp_class());
	

	std::chrono::high_resolution_clock::time_point vpd_input_0 = std::chrono::high_resolution_clock::now();
	bool input_0_verify = input_vpd::verify(r_0_mpz, digest_input.first[0], digest_input.second[0], p -> Zu.to_gmp_class(), input_0_mpz, witnesses_0.first, witnesses_0.second);
	std::chrono::high_resolution_clock::time_point vpd_input_1 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> vpd_input_span = std::chrono::duration_cast<std::chrono::duration<double>>(vpd_input_1 - vpd_input_0);
	verification_time += vpd_input_span.count();

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
		std::cerr << "Verification Time " << verification_time << std::endl;
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
	for(int i = 0; i < C.total_depth; ++i)
	{
		delete[] C.circuit[i].gates;
	}
	delete[] C.circuit;
}
