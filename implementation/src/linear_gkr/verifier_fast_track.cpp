#include "linear_gkr/verifier_fast_track.h"
#include <string>
#include <utility>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include "linear_gkr/random_generator.h"

void verifier::get_prover(prover *pp)
{
	p = pp;
}

void verifier::read_circuit(const char *path)
{
	int d;
	static char str[300];
	FILE *circuit_in;
	circuit_in = fopen(path, "r");

	fscanf(circuit_in, "%d", &d);
	int n;
	C.circuit = new layer[d];
	C.total_depth = d;
	int max_bit_length = -1;
	for(int i = 0; i < d; ++i)
	{
		fscanf(circuit_in, "%d", &n);
		if(n == 1)
			C.circuit[i].gates = new gate[2];
		else
			C.circuit[i].gates = new gate[n];
		int max_gate = -1;
		int previous_g = -1;
		for(int j = 0; j < n; ++j)
		{
			int ty, g, u, v;
			fscanf(circuit_in, "%d%d%d%d", &ty, &g, &u, &v);
			if(g != previous_g + 1)
			{
				printf("Error, gates must be in sorted order, and full [0, 2^n - 1].");
			}
			previous_g = g;
			C.circuit[i].gates[g] = gate(ty, u, v);
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
			C.circuit[i].gates[max_gate] = gate(2, 0, 0);
			C.circuit[i].bit_length = cnt;
		}
		else
		{
			C.circuit[i].bit_length = cnt - 1;
		}
		fprintf(stderr, "layer %d, bit_length %d\n", i, C.circuit[i].bit_length);
		if(C.circuit[i].bit_length > max_bit_length)
			max_bit_length = C.circuit[i].bit_length;
	}
	p -> init_array(max_bit_length);
	fclose(circuit_in);
}

prime_field::field_element verifier::add(int depth, 
	const prime_field::field_element* z, const prime_field::field_element* r_u, const prime_field::field_element* r_v, const prime_field::field_element* one_minus_z, const prime_field::field_element* one_minus_r_u, const prime_field::field_element* one_minus_r_v)
{
	//brute force for sanity check
	//it's slow
	prime_field::field_element ret = prime_field::field_element(0);
	for(int i = 0; i < (1 << C.circuit[depth].bit_length); ++i)
	{
		int g = i;
		prime_field::field_element cur = prime_field::field_element(1);
		if(C.circuit[depth].gates[g].ty == 0)
		{
			int u, v;
			u = C.circuit[depth].gates[g].u;
			v = C.circuit[depth].gates[g].v;
			for(int j = 0; j < C.circuit[depth].bit_length; ++j)
			{
				if((g & 1) == 0)
				{
					cur.value = cur.value * one_minus_z[j].value;
				}
				else
				{
					cur.value = cur.value * z[j].value;
				}
				g >>= 1;
			}
			for(int j = 0; j < C.circuit[depth - 1].bit_length; ++j)
			{
				if((u & 1) == 0)
				{
					cur.value = cur.value * one_minus_r_u[j].value;
				}
				else
				{
					cur.value = cur.value * r_u[j].value;
				}
				u >>= 1;
			}
			for(int j = 0; j < C.circuit[depth - 1].bit_length; ++j)
			{
				if((v & 1) == 0)
				{
					cur.value = cur.value * one_minus_r_v[j].value;
				}
				else
				{
					cur.value = cur.value * r_v[j].value;
				}
				v >>= 1;
			}
			ret.value = ret.value + cur.value;
		}
	}
	ret.value = ret.value % prime_field::mod;
	return ret;
}
prime_field::field_element verifier::mult(int depth, 
	const prime_field::field_element* z, const prime_field::field_element* r_u, const prime_field::field_element* r_v, const prime_field::field_element* one_minus_z, const prime_field::field_element* one_minus_r_u, const prime_field::field_element* one_minus_r_v)
{
	//also brute force
	prime_field::field_element ret = prime_field::field_element(0);
	for(int i = 0; i < (1 << C.circuit[depth].bit_length); ++i)
	{
		int g = i;
		prime_field::field_element cur = prime_field::field_element(1);
		if(C.circuit[depth].gates[g].ty == 1)
		{
			int u, v;
			u = C.circuit[depth].gates[g].u;
			v = C.circuit[depth].gates[g].v;
			for(int j = 0; j < C.circuit[depth].bit_length; ++j)
			{
				if((g & 1) == 0)
				{
					cur.value = cur.value * one_minus_z[j].value;
				}
				else
				{
					cur.value = cur.value * z[j].value;
				}
				g >>= 1;
			}
			for(int j = 0; j < C.circuit[depth - 1].bit_length; ++j)
			{
				if((u & 1) == 0)
				{
					cur.value = cur.value * one_minus_r_u[j].value;
				}
				else
				{
					cur.value = cur.value * r_u[j].value;
				}
				u >>= 1;
			}
			for(int j = 0; j < C.circuit[depth - 1].bit_length; ++j)
			{
				if((v & 1) == 0)
				{
					cur.value = cur.value * one_minus_r_v[j].value;
				}
				else
				{
					cur.value = cur.value * r_v[j].value;
				}
				v >>= 1;
			}
			ret.value = ret.value + cur.value;
		}
	}
	ret.value = ret.value % prime_field::mod;
	return ret;
}

gmp_randstate_t rstate;

prime_field::field_element* generate_randomness(unsigned int size)
{
	int k = size;
	prime_field::field_element* ret;
	ret = new prime_field::field_element[k];

	for(int i = 0; i < k; ++i)
	{
		mpz_t random_element;
		mpz_init(random_element);
		mpz_urandomm(random_element, rstate, prime_field::mod.get_mpz_t());
		ret[i] = prime_field::field_element(mpz_class(random_element));
	}
	return ret;
}

prime_field::field_element verifier::V_in(const prime_field::field_element* r_0, const prime_field::field_element* one_minus_r_0,
								prime_field::field_element* output_raw, int r_0_size, int output_size)
{
	prime_field::field_element* output = new prime_field::field_element[output_size];
	for(int i = 0; i < output_size; ++i)
		output[i] = output_raw[i];
	for(int i = 0; i < r_0_size; ++i)
	{
		int last_gate;
		int cnt = 0;
		for(int j = 0; j < (output_size >> 1); ++j)
			output[j] = output[j << 1] * (one_minus_r_0[i]) + output[j << 1 | 1] * (r_0[i]);
	}
	auto ret = output[0];
	delete[] output;
	return ret;
}

bool verifier::verify()
{
	gmp_randinit_default(rstate);
	p -> proof_init();

	auto result = p -> evaluate();
	fprintf(stderr, "evaluation result:\n");
//	for(int i = 0; i < (1 << C.circuit[C.total_depth - 1].bit_length); ++i)
//	{
//		fprintf(stderr, "%d %s\n", i, result[i].to_string(10).c_str());
//	}
//	fprintf(stderr, "\n");

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
	prime_field::field_element a_0 = p -> V_res(one_minus_r_0, r_0, result, C.circuit[C.total_depth - 1].bit_length, (1 << (C.circuit[C.total_depth - 1].bit_length)));
	a_0 = alpha * a_0;
	prime_field::field_element a_1 = prime_field::field_element(0); //* beta

	printf("a_0 = %s\n", a_0.to_string(10).c_str());

	prime_field::field_element alpha_beta_sum = a_0; //+ a_1

	for(int i = C.total_depth - 1; i >= 1; --i)
	{
		p -> sumcheck_init(i, C.circuit[i].bit_length, C.circuit[i - 1].bit_length, C.circuit[i - 1].bit_length, alpha, beta, r_0, r_1, one_minus_r_0, one_minus_r_1);
		p -> sumcheck_phase1_init();
		prime_field::field_element previous_random = prime_field::field_element(0);
		//next level random
		auto r_u = generate_randomness(C.circuit[i - 1].bit_length);
		auto r_v = generate_randomness(C.circuit[i - 1].bit_length);
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
			quadratic_poly poly = p -> sumcheck_phase1_update(previous_random, j);
			previous_random = r_u[j];
			if(poly.eval(0) + poly.eval(1) != alpha_beta_sum)
			{
				fprintf(stderr, "Verification fail, phase1, circuit %d, current bit %d\n", i, j);
				return false;
			}
			else
			{
				fprintf(stderr, "Verification Pass, phase1, circuit %d, current bit %d\n", i, j);
			}
			alpha_beta_sum = poly.eval(r_u[j]);
		}
//		std::cout << "phase1 sum " << alpha_beta_sum.to_string(10) << std::endl;
		p -> sumcheck_phase2_init(previous_random, r_u, one_minus_r_u);
		previous_random = prime_field::field_element(0);
		for(int j = 0; j < C.circuit[i - 1].bit_length; ++j)
		{
			quadratic_poly poly = p -> sumcheck_phase2_update(previous_random, j);
			previous_random = r_v[j];
			if(poly.eval(0) + poly.eval(1) != alpha_beta_sum)
			{
				fprintf(stderr, "Verification fail, phase2, circuit level %d, current bit %d\n", i, j);
				return false;
			}
			else
			{
				fprintf(stderr, "Verification Pass, phase2, circuit level %d, current bit %d\n", i, j);
			}
			alpha_beta_sum = poly.eval(r_v[j]);
		}
//		std::cout << "phase2 sum " << alpha_beta_sum.to_string(10) << std::endl;
		auto final_claims = p -> sumcheck_finalize(previous_random);
		auto v_u = final_claims.first;
		auto v_v = final_claims.second;

//		std::cout << "v_u = " << v_u.to_string(10) << std::endl;
//		std::cout << "v_v = " << v_v.to_string(10) << std::endl;

//		std::cout << "alpha = " << alpha.to_string(10) << std::endl;
//		std::cout << "beta = " << beta.to_string(10) << std::endl;

		auto mult_value = mult(i, r_0, r_u, r_v, one_minus_r_0, one_minus_r_u, one_minus_r_v) * alpha + 
						mult(i, r_1, r_u, r_v, one_minus_r_1, one_minus_r_u, one_minus_r_v) * beta;
		auto add_value = add(i, r_0, r_u, r_v, one_minus_r_0, one_minus_r_u, one_minus_r_v) * alpha + 
						add(i, r_1, r_u, r_v, one_minus_r_1, one_minus_r_u, one_minus_r_v) * beta;
//		std::cout << "mult_value = " << mult_value.to_string(10) << std::endl;
//		std::cout << "add_value = " << add_value.to_string(10) << std::endl;

		if(alpha_beta_sum != add_value * (v_u + v_v) + mult_value * v_u * v_v)
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
		alpha_beta_sum = alpha * v_u + beta * v_v;

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
	prime_field::field_element* input;
	input = new prime_field::field_element[(1 << C.circuit[0].bit_length)];

	for(int i = 0; i < (1 << C.circuit[0].bit_length); ++i)
	{
		int g = i;
		if(C.circuit[0].gates[g].ty == 3)
		{
			input[g] = prime_field::field_element(C.circuit[0].gates[g].u);
		}
		else
			assert(false);
	}
	auto input_0 = V_in(r_0, one_minus_r_0, input, C.circuit[0].bit_length, (1 << C.circuit[0].bit_length)), 
		 input_1 = V_in(r_1, one_minus_r_1, input, C.circuit[0].bit_length, (1 << C.circuit[0].bit_length));

	delete[] input;
	delete[] r_0;
	delete[] r_1;
	delete[] one_minus_r_0;
	delete[] one_minus_r_1;
	if(alpha_beta_sum != input_0 * alpha + input_1 * beta)
	{
		fprintf(stderr, "Verification fail, final input check fail.\n");
		return false;
	}
	else
	{
		fprintf(stderr, "Verification pass\n");
		fprintf(stderr, "Prove Time %f\n", (float)(p -> total_time) / (float)CLOCKS_PER_SEC);
	}
	p -> delete_self();
	delete_self();
	return true;
}

void verifier::delete_self()
{
	for(int i = 0; i < C.total_depth; ++i)
	{
		delete[] C.circuit[i].gates;
	}
	delete[] C.circuit;
}