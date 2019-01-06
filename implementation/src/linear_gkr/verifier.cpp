#include "linear_gkr/verifier.h"
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
	C.circuit.clear();
	for(int i = 0; i < d; ++i)
	{
		C.circuit.push_back(layer());
		fscanf(circuit_in, "%d", &n);
		int max_gate = -1;
		int previous_g = -1;
		for(int j = 0; j < n; ++j)
		{
			int ty, g, u, v;
			fscanf(circuit_in, "%d%d%d%d", &ty, &g, &u, &v);
			if(g < previous_g)
			{
				printf("Error, gates must be in sorted order.");
			}
			assert(g > previous_g);
			previous_g = g;
			if(g > max_gate)
				max_gate = g;
			C.circuit[i].gates[g] = std::make_pair(ty, std::make_pair(u, v));
			C.circuit[i].gate_id.push_back(g);
		}
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
			C.circuit[i].gate_id.push_back(max_gate);
			C.circuit[i].gates[max_gate] = std::make_pair(2, std::make_pair(0, 0));
			C.circuit[i].bit_length = cnt;
		}
		else
		{
			C.circuit[i].bit_length = cnt - 1;
		}
		fprintf(stderr, "layer %d, bit_length %d\n", i, C.circuit[i].bit_length);
	}

	fclose(circuit_in);
}

prime_field::field_element verifier::add(int depth, 
	std::vector<prime_field::field_element> z, std::vector<prime_field::field_element> r_u, std::vector<prime_field::field_element> r_v)
{
	//brute force for sanity check
	//it's slow
	prime_field::field_element ret = prime_field::field_element(0);
	for(int i = 0; i < C.circuit[depth].gate_id.size(); ++i)
	{
		int g = C.circuit[depth].gate_id[i];
		prime_field::field_element cur = prime_field::field_element(1);
		if(C.circuit[depth].gates[g].first == 0)
		{
			int u, v;
			u = C.circuit[depth].gates[g].second.first;
			v = C.circuit[depth].gates[g].second.second;
			for(int j = 0; j < C.circuit[depth].bit_length; ++j)
			{
				if((g & 1) == 0)
				{
					cur = cur * (prime_field::field_element(1) - z[j]);
				}
				else
				{
					cur = cur * z[j];
				}
				g >>= 1;
			}
			for(int j = 0; j < C.circuit[depth - 1].bit_length; ++j)
			{
				if((u & 1) == 0)
				{
					cur = cur * (prime_field::field_element(1) - r_u[j]);
				}
				else
				{
					cur = cur * r_u[j];
				}
				u >>= 1;
			}
			for(int j = 0; j < C.circuit[depth - 1].bit_length; ++j)
			{
				if((v & 1) == 0)
				{
					cur = cur * (prime_field::field_element(1) - r_v[j]);
				}
				else
				{
					cur = cur * r_v[j];
				}
				v >>= 1;
			}
			ret = ret + cur;
		}
	}
	return ret;
}
prime_field::field_element verifier::mult(int depth, 
	std::vector<prime_field::field_element> z, std::vector<prime_field::field_element> r_u, std::vector<prime_field::field_element> r_v)
{
	//also brute force
	prime_field::field_element ret = prime_field::field_element(0);
	for(int i = 0; i < C.circuit[depth].gate_id.size(); ++i)
	{
		int g = C.circuit[depth].gate_id[i];
		prime_field::field_element cur = prime_field::field_element(1);
		if(C.circuit[depth].gates[g].first == 1)
		{
			int u, v;
			u = C.circuit[depth].gates[g].second.first;
			v = C.circuit[depth].gates[g].second.second;
			for(int j = 0; j < C.circuit[depth].bit_length; ++j)
			{
				if((g & 1) == 0)
				{
					cur = cur * (prime_field::field_element(1) - z[j]);
				}
				else
				{
					cur = cur * z[j];
				}
				g >>= 1;
			}
			for(int j = 0; j < C.circuit[depth - 1].bit_length; ++j)
			{
				if((u & 1) == 0)
				{
					cur = cur * (prime_field::field_element(1) - r_u[j]);
				}
				else
				{
					cur = cur * r_u[j];
				}
				u >>= 1;
			}
			for(int j = 0; j < C.circuit[depth - 1].bit_length; ++j)
			{
				if((v & 1) == 0)
				{
					cur = cur * (prime_field::field_element(1) - r_v[j]);
				}
				else
				{
					cur = cur * r_v[j];
				}
				v >>= 1;
			}
			ret = ret + cur;
		}
	}
	return ret;
}

gmp_randstate_t rstate;

std::vector<prime_field::field_element> generate_randomness(unsigned int size)
{
	int k = size;
	std::vector<prime_field::field_element> ret;

	for(int i = 0; i < k; ++i)
	{
		mpz_t random_element;
		mpz_init(random_element);
		mpz_urandomm(random_element, rstate, prime_field::mod.get_mpz_t());
		ret.push_back(mpz_class(random_element));
	}
	return ret;
}

prime_field::field_element verifier::V_in(const std::vector<prime_field::field_element> &r_0, 
								std::vector<std::pair<int, prime_field::field_element> > output)
{
	std::vector<std::pair<int, prime_field::field_element> > tmp;
	for(int i = 0; i < r_0.size(); ++i)
	{
		tmp.clear();
		int last_gate;
		int cnt = 0;
		for(int j = 0; j < output.size(); ++j)
		{
			prime_field::field_element m = r_0[i];
			if((output[j].first & 1) == 0)
			{
				m = (prime_field::field_element(mpz_class(1)) - m);
			}
			if(j == 0)
			{
				tmp.push_back(std::make_pair(output[j].first >> 1, output[j].second * m));
				last_gate = output[j].first >> 1;
				cnt++;
			}
			else
			{
				if((output[j].first >> 1) == last_gate)
				{
					tmp[cnt - 1] = std::make_pair(last_gate, tmp[cnt - 1].second + output[j].second * m);
				}
				else
				{
					last_gate = output[j].first >> 1;
					tmp.push_back(std::make_pair(output[j].first >> 1, output[j].second * m));
					cnt++;
				}
			}
		}
		output = tmp;
	}
	assert(output.size() == 1);
	return output[0].second;
}

bool verifier::verify()
{
	gmp_randinit_default(rstate);
	p -> proof_init();

	auto result = p -> evaluate();
//	fprintf(stderr, "evaluation result:\n");
//	for(auto x : result)
//	{
//		fprintf(stderr, "%d %s\n", x.first, x.second.to_string(10).c_str());
//	}
//	fprintf(stderr, "\n");

	prime_field::field_element alpha, beta;
	alpha.value = 1;
	beta.value = 0;
	random_oracle oracle;
	//initial random value
	std::vector<prime_field::field_element> r_0 = generate_randomness(C.circuit[C.circuit.size() - 1].bit_length), r_1 = generate_randomness(C.circuit[C.circuit.size() - 1].bit_length);

	auto a_0 = p -> V_res(r_0, result);
	a_0 = alpha * a_0;
	prime_field::field_element a_1 = prime_field::field_element(mpz_class(0)) * beta;

	printf("a_0 = %s\n", a_0.to_string(10).c_str());

	auto alpha_beta_sum = a_0 + a_1;

	for(int i = (int)C.circuit.size() - 1; i >= 1; --i)
	{
		p -> sumcheck_init(i, C.circuit[i].bit_length, C.circuit[i - 1].bit_length, C.circuit[i - 1].bit_length, alpha, beta, r_0, r_1);
		p -> sumcheck_phase1_init();
		prime_field::field_element previous_random = prime_field::field_element(0);
		//next level random
		auto r_u = generate_randomness(C.circuit[i - 1].bit_length);
		auto r_v = generate_randomness(C.circuit[i - 1].bit_length);

//		for(int i = 0; i < r_u.size(); ++i)
//			std::cout << "r_u[" << i << "] = " << r_u[i].to_string(10) << std::endl;

//		for(int i = 0; i < r_v.size(); ++i)
//			std::cout << "r_v[" << i << "] = " << r_v[i].to_string(10) << std::endl;


//		for(int i = 0; i < r_0.size(); ++i)
//			std::cout << "r_0[" << i << "] = " << r_0[i].to_string(10) << std::endl;

//		for(int i = 0; i < r_1.size(); ++i)
//			std::cout << "r_1[" << i << "] = " << r_1[i].to_string(10) << std::endl;


		for(int j = 0; j < C.circuit[i - 1].bit_length; ++j)
		{
			quadratic_poly poly = p -> sumcheck_phase1_update(previous_random);
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
		std::cout << "phase1 sum " << alpha_beta_sum.to_string(10) << std::endl;
		p -> sumcheck_phase2_init(previous_random, r_u);
		previous_random = prime_field::field_element(0);
		for(int j = 0; j < C.circuit[i - 1].bit_length; ++j)
		{
			quadratic_poly poly = p -> sumcheck_phase2_update(previous_random);
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
		std::cout << "phase2 sum " << alpha_beta_sum.to_string(10) << std::endl;
		auto final_claims = p -> sumcheck_finalize(previous_random);
		auto v_u = final_claims.first;
		auto v_v = final_claims.second;

		std::cout << "v_u = " << v_u.to_string(10) << std::endl;
		std::cout << "v_v = " << v_v.to_string(10) << std::endl;

		std::cout << "alpha = " << alpha.to_string(10) << std::endl;
		std::cout << "beta = " << beta.to_string(10) << std::endl;

		auto mult_value = mult(i, r_0, r_u, r_v) * alpha + mult(i, r_1, r_u, r_v) * beta;
		auto add_value = add(i, r_0, r_u, r_v) * alpha + add(i, r_1, r_u, r_v) * beta;
		std::cout << "mult_value = " << mult_value.to_string(10) << std::endl;
		std::cout << "add_value = " << add_value.to_string(10) << std::endl;

		if(alpha_beta_sum != add_value * (v_u + v_v) + mult_value * v_u * v_v)
		{
			fprintf(stderr, "Verification fail, semi final, circuit level %d\n", i);
			return false;
		}
		else
		{
			fprintf(stderr, "Verification Pass, semi final, circuit level %d\n", i);
		}
		alpha = generate_randomness(1)[0];
		beta = generate_randomness(1)[0];
		alpha_beta_sum = alpha * v_u + beta * v_v;
		r_0 = r_u;
		r_1 = r_v;
		//todo randomize alpha beta
	}

	//post sumcheck
	std::vector<std::pair<int, prime_field::field_element> > input;
	for(int i = 0; i < C.circuit[0].gate_id.size(); ++i)
	{
		int g = C.circuit[0].gate_id[i];
		if(C.circuit[0].gates[g].first == 3)
		{
			input.push_back(std::make_pair(g, prime_field::field_element(C.circuit[0].gates[g].second.first)));
		}
	}
	auto input_0 = V_in(r_0, input), input_1 = V_in(r_1, input);
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
	return true;
}