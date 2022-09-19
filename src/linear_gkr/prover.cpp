#include "linear_gkr/prover.h"

void prover::get_circuit(const layered_circuit &from_verifier)
{
	C = from_verifier;
}

prime_field::field_element prover::V_res(const std::vector<prime_field::field_element> &r_0, 
								std::vector<std::pair<int, prime_field::field_element> > output)
{
	clock_t t0 = clock();
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
	total_time += (clock() - t0);
	return output[0].second;
}

std::vector<std::pair<int, prime_field::field_element> > prover::evaluate()
{
	clock_t t0 = clock();
	circuit_value.clear();
	circuit_value.push_back(std::unordered_map<int, prime_field::field_element >());
	for(int i = 0; i < C.circuit[0].gate_id.size(); ++i)
	{
		int g, u, v, ty;
		g = C.circuit[0].gate_id[i];
		std::pair<int, std::pair<int, int> > info = C.circuit[0].gates[g];
		ty = info.first;
		u = info.second.first;
		v = info.second.second;
		assert(ty == 3);
		circuit_value[0][g] = mpz_class(u);
	}
	std::vector<std::pair<int, prime_field::field_element> > ret;
	for(int i = 1; i < C.circuit.size(); ++i)
	{
		circuit_value.push_back(std::unordered_map<int, prime_field::field_element>());
		for(int j = 0; j < C.circuit[i].gate_id.size(); ++j)
		{
			int g, u, v, ty;
			g = C.circuit[i].gate_id[j];
			std::pair<int, std::pair<int, int> > info = C.circuit[i].gates[g];
			ty = info.first;
			u = info.second.first;
			v = info.second.second;
			if(ty == 0)
			{
				circuit_value[i][g] = circuit_value[i - 1][u] + circuit_value[i - 1][v];
			}
			else if(ty == 1)
			{
				circuit_value[i][g] = circuit_value[i - 1][u] * circuit_value[i - 1][v];
			}
			else if(ty == 2)
			{
				circuit_value[i][g] = mpz_class(0);
			}
			else if(ty == 3)
			{
				circuit_value[i][g] = mpz_class(u);
			}
			else
			{
				assert(false);
			}
			if(i + 1 == C.circuit.size())
			{
				ret.push_back(std::make_pair(g, circuit_value[i][g]));
			}
		}
	}
/*
	//cheating prover test
	ret[0].second =  prime_field::field_element(111);
*/
	fprintf(stderr, "total evaluation time: %f\n", ((float)clock() - (float)t0) / CLOCKS_PER_SEC);
	return ret;
}

void prover::sumcheck_init(int layer_id, int bit_length_g, int bit_length_u, int bit_length_v, 
	const prime_field::field_element &a, const prime_field::field_element &b, const std::vector<prime_field::field_element> &R_0,
	const std::vector<prime_field::field_element> &R_1)
{
	r_0 = R_0;
	r_1 = R_1;
	alpha = a;
	beta = b;
	randomness_from_verifier.clear();
	sumcheck_layer_id = layer_id;
	length_g = bit_length_g;
	length_u = bit_length_u;
	length_v = bit_length_v;
}
void prover::DFS(std::unordered_map<int, linear_poly> &A, int g, int depth, 
	prime_field::field_element alpha_value, prime_field::field_element beta_value, const int gate_type)
{
	if(depth == length_g)
	{
		int u, v;
		if(C.circuit[sumcheck_layer_id].gates[g].first != gate_type) //not a mult gate
		{
			return;
		}
		u = C.circuit[sumcheck_layer_id].gates[g].second.first;
		v = C.circuit[sumcheck_layer_id].gates[g].second.second;
		A[u] = A[u] + circuit_value[sumcheck_layer_id - 1][v] * (alpha_value * alpha + beta_value * beta);
	}
	else
	{
		g &= ((-1) ^ (1 << depth));
		DFS(A, g, depth + 1, 
			alpha_value * (prime_field::field_element(1) - r_0[depth]), beta_value * (prime_field::field_element(1) - r_1[depth]), gate_type);
		g |= (1 << depth);
		DFS(A, g, depth + 1, alpha_value * r_0[depth], beta_value * r_1[depth], gate_type);
	}
}

void prover::DFS_add(std::unordered_map<int, linear_poly> &A, int g, int depth,
					prime_field::field_element alpha_value, prime_field::field_element beta_value, const int gate_type)
{
	if(depth == length_g)
	{
		int u, v;
		if(C.circuit[sumcheck_layer_id].gates[g].first != gate_type) //not a mult gate
		{
			return;
		}
		u = C.circuit[sumcheck_layer_id].gates[g].second.first;
		v = C.circuit[sumcheck_layer_id].gates[g].second.second;
		A[u] = A[u] + alpha_value * alpha + beta_value * beta;
	}
	else
	{
		g &= ((-1) ^ (1 << depth));
		DFS_add(A, g, depth + 1, 
			alpha_value * (prime_field::field_element(1) - r_0[depth]), beta_value * (prime_field::field_element(1) - r_1[depth]), gate_type);
		g |= (1 << depth);
		DFS_add(A, g, depth + 1, alpha_value * r_0[depth], beta_value * r_1[depth], gate_type);
	}
}
void prover::sumcheck_phase1_init()
{
	fprintf(stderr, "sumcheck level %d, phase1 init start\n", sumcheck_layer_id);
	randomness_from_verifier.clear();
	clock_t t0 = clock();
	//mult init
	mult_array.clear();
	V_mult.clear();
	for(int i = 0; i < C.circuit[sumcheck_layer_id - 1].gate_id.size(); ++i)
	{
		int g = C.circuit[sumcheck_layer_id - 1].gate_id[i];
		mult_array[g] = 
			linear_poly(prime_field::field_element(0), prime_field::field_element(0));
		V_mult[g] = 
			circuit_value[sumcheck_layer_id - 1][g];
	}
	DFS(mult_array, 0, 0, prime_field::field_element(1), prime_field::field_element(1), 1);
	//add init
	addV_array.clear();
	add_array.clear();
	V_add.clear();

	for(int i = 0; i < C.circuit[sumcheck_layer_id - 1].gate_id.size(); ++i)
	{
		int g = C.circuit[sumcheck_layer_id - 1].gate_id[i];
		addV_array[g] = 
			linear_poly(prime_field::field_element(0), prime_field::field_element(0));
		add_array[g] = 
			linear_poly(prime_field::field_element(0), prime_field::field_element(0));
		V_add[C.circuit[sumcheck_layer_id - 1].gate_id[i]] = 
			circuit_value[sumcheck_layer_id - 1][g];
	}
	DFS(addV_array, 0, 0, prime_field::field_element(1), prime_field::field_element(1), 0);
	DFS_add(add_array, 0, 0, prime_field::field_element(1), prime_field::field_element(1), 0);
	current_sumcheck_gates = C.circuit[sumcheck_layer_id - 1].gate_id;
	fprintf(stderr, "sumcheck level %d, phase1 init finished\n", sumcheck_layer_id);
	total_time += clock() - t0;
}

quadratic_poly prover::sumcheck_phase1_update(prime_field::field_element previous_random)
{
	clock_t t0 = clock();
	std::vector<int> nxt_gate_id;
	for(int i = 0; i < current_sumcheck_gates.size(); ++i)
	{
		int g = current_sumcheck_gates[i];
		int target_g = (g >> 1);
		nxt_gate_id.push_back(target_g);
	}
	nxt_gate_id.resize(std::distance(nxt_gate_id.begin(), std::unique(nxt_gate_id.begin(), nxt_gate_id.end())));
	quadratic_poly ret = quadratic_poly(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));
	std::unordered_map<int, linear_poly> mult_addV, V, new_add;
	for(int i = 0; i < nxt_gate_id.size(); ++i)
	{
		prime_field::field_element zero_value, one_value;
		int g = nxt_gate_id[i];
		int g_zero = g << 1, g_one = g << 1 | 1;
		if(mult_array.find(g_zero) != mult_array.end())
		{
			zero_value = mult_array[g_zero].eval(previous_random);
		}
		else
		{
			zero_value = prime_field::field_element(0);
		}
		if(mult_array.find(g_one) != mult_array.end())
		{
			one_value = mult_array[g_one].eval(previous_random);
		}
		else
		{
			one_value = prime_field::field_element(0);
		}
		mult_addV[g] = linear_poly(one_value - zero_value, zero_value);

		if(V_mult.find(g_zero) != V_mult.end())
		{
			zero_value = V_mult[g_zero].eval(previous_random);
		}
		else
		{
			zero_value = prime_field::field_element(0);
		}

		if(V_mult.find(g_one) != V_mult.end())
		{
			one_value = V_mult[g_one].eval(previous_random);
		}
		else
		{
			one_value = prime_field::field_element(0);
		}

		V[g] = linear_poly(one_value - zero_value, zero_value);
		ret = ret + mult_addV[g] * V[g];
	}
	//add gate

	mult_array = mult_addV;
	V_mult = V;
	mult_addV.clear(), V.clear();
	new_add.clear();
	for(int i = 0; i < nxt_gate_id.size(); ++i)
	{
		prime_field::field_element zero_value, one_value;
		int g = nxt_gate_id[i];
		int g_zero = g << 1, g_one = g << 1 | 1;
		if(addV_array.find(g_zero) != addV_array.end())
		{
			zero_value = addV_array[g_zero].eval(previous_random);
		}
		else
		{
			zero_value = prime_field::field_element(0);
		}
		if(addV_array.find(g_one) != addV_array.end())
		{
			one_value = addV_array[g_one].eval(previous_random);
		}
		else
		{
			one_value = prime_field::field_element(0);
		}
		mult_addV[g] = linear_poly(one_value - zero_value, zero_value);
		if(V_add.find(g_zero) != V_add.end())
		{
			zero_value = V_add[g_zero].eval(previous_random);
		}
		else
		{
			zero_value = prime_field::field_element(0);
		}
		if(V_add.find(g_one) != V_add.end())
		{
			one_value = V_add[g_one].eval(previous_random);
		}
		else
		{
			one_value = prime_field::field_element(0);
		}
		V[g] = linear_poly(one_value - zero_value, zero_value);

		if(add_array.find(g_zero) != add_array.end())
		{
			zero_value = add_array[g_zero].eval(previous_random);
		}
		else
		{
			zero_value = prime_field::field_element(0);
		}
		if(add_array.find(g_one) != add_array.end())
		{
			one_value = add_array[g_one].eval(previous_random);
		}
		else
		{
			one_value = prime_field::field_element(0);
		}
		new_add[g] = linear_poly(one_value - zero_value, zero_value);
		ret = ret + new_add[g] * V[g] + quadratic_poly(0, mult_addV[g].a, mult_addV[g].b);
	}
	addV_array = mult_addV;
	add_array = new_add;
	V_add = V;
	current_sumcheck_gates = nxt_gate_id;
	total_time += clock() - t0;
	return ret;
}

void prover::DFS_betag(std::unordered_map<int, prime_field::field_element> &A, std::vector<prime_field::field_element> r,
	int g, int depth, prime_field::field_element value, int gate_type)
{
	if(depth == length_g)
	{
		A[g] = value;
	}
	else
	{
		g &= ((-1) ^ (1 << depth));
		DFS_betag(A, r, g, depth + 1, 
			value * (prime_field::field_element(1) - r[depth]), gate_type);
		g |= (1 << depth);
		DFS_betag(A, r, g, depth + 1, value * r[depth], gate_type);
	}
}

void prover::DFS_betau(std::unordered_map<int, prime_field::field_element> &A, std::vector<prime_field::field_element> r,
	int u, int depth, prime_field::field_element value, int gate_type)
{
	if(depth == length_u)
	{
		A[u] = value;
	}
	else
	{
		u &= ((-1) ^ (1 << depth));
		DFS_betau(A, r, u, depth + 1, value * (prime_field::field_element(1) - r[depth]), gate_type);
		u |= (1 << depth);
		DFS_betau(A, r, u, depth + 1, value * r[depth], gate_type);
	}
}

void prover::sumcheck_phase2_init(prime_field::field_element previous_random, std::vector<prime_field::field_element> r_u)
{
	fprintf(stderr, "sumcheck level %d, phase2 init start\n", sumcheck_layer_id);
	clock_t t0 = clock();
	v_u = (*V_add.begin()).second.eval(previous_random);
//	fprintf(stderr, "v_u %s\n", v_u.to_string(10).c_str());
	randomness_from_verifier.clear();
	
	//mult
	//init betag
	std::unordered_map<int, prime_field::field_element> beta_g_r0, beta_g_r1;
	DFS_betag(beta_g_r0, r_0, 0, 0, prime_field::field_element(1), 1);
	DFS_betag(beta_g_r1, r_1, 0, 0, prime_field::field_element(1), 1);
	//init betau
	std::unordered_map<int, prime_field::field_element> beta_u;
	DFS_betau(beta_u, r_u, 0, 0, prime_field::field_element(1), 1);

	mult_array.clear();
	V_mult.clear();
	add_array.clear();
	addV_array.clear();
	V_add.clear();

	for(int i = 0; i < C.circuit[sumcheck_layer_id - 1].gate_id.size(); ++i)
	{
		int v = C.circuit[sumcheck_layer_id - 1].gate_id[i];
		mult_array[v] = linear_poly(prime_field::field_element(0), prime_field::field_element(0));
		add_array[v] = linear_poly(prime_field::field_element(0), prime_field::field_element(0));
		addV_array[v] = linear_poly(prime_field::field_element(0), prime_field::field_element(0));
		V_mult[v] = circuit_value[sumcheck_layer_id - 1][v];
		V_add[v] = circuit_value[sumcheck_layer_id - 1][v];
	}

	for(int i = 0; i < C.circuit[sumcheck_layer_id].gate_id.size(); ++i)
	{
		int g = C.circuit[sumcheck_layer_id].gate_id[i];
		int ty = C.circuit[sumcheck_layer_id].gates[g].first;
		int u = C.circuit[sumcheck_layer_id].gates[g].second.first;
		int v = C.circuit[sumcheck_layer_id].gates[g].second.second;
		if(ty == 1) //mult gate
		{
			mult_array[v] = mult_array[v] + linear_poly((beta_g_r0[g] * beta_u[u] * alpha + beta_g_r1[g] * beta_u[u] * beta) * v_u);
		}
	}

	beta_g_r0.clear(), beta_g_r1.clear();
	beta_u.clear();

	DFS_betag(beta_g_r0, r_0, 0, 0, prime_field::field_element(1), 0);
	DFS_betag(beta_g_r1, r_1, 0, 0, prime_field::field_element(1), 0);
	DFS_betau(beta_u, r_u, 0, 0, prime_field::field_element(1), 0);

	for(int i = 0; i < C.circuit[sumcheck_layer_id].gate_id.size(); ++i)
	{
		int g = C.circuit[sumcheck_layer_id].gate_id[i];
		int ty = C.circuit[sumcheck_layer_id].gates[g].first;
		int u = C.circuit[sumcheck_layer_id].gates[g].second.first;
		int v = C.circuit[sumcheck_layer_id].gates[g].second.second;
		if(ty == 0) //add gate
		{
			add_array[v] = add_array[v] + linear_poly(beta_g_r0[g] * beta_u[u] * alpha + beta_g_r1[g] * beta_u[u] * beta);
			addV_array[v] = addV_array[v] + linear_poly((beta_g_r0[g] * beta_u[u] * alpha + beta_g_r1[g] * beta_u[u] * beta) * v_u);
		}
	}
//	for(int i = 0; i < (1 << C.circuit[sumcheck_layer_id - 1].bit_length); ++i)
//	{
//		fprintf(stderr, "add %s\n", add_array[i].b.to_string(10).c_str());
//		fprintf(stderr, "addV %s\n", add_array[i].b.to_string(10).c_str());
//		fprintf(stderr, "V_mult_add %s\n", add_array[i].b.to_string(10).c_str());
//	}
	current_sumcheck_gates = C.circuit[sumcheck_layer_id - 1].gate_id;
	total_time += clock() - t0;
}

quadratic_poly prover::sumcheck_phase2_update(prime_field::field_element previous_random)
{
	clock_t t0 = clock();
	std::vector<int> nxt_gate_id;
	for(int i = 0; i < current_sumcheck_gates.size(); ++i)
	{
		int g = current_sumcheck_gates[i];
		int target_g = (g >> 1);
		nxt_gate_id.push_back(target_g);
	}
	nxt_gate_id.resize(std::distance(nxt_gate_id.begin(), std::unique(nxt_gate_id.begin(), nxt_gate_id.end())));
	quadratic_poly ret = quadratic_poly(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));
	std::unordered_map<int, linear_poly> mult_addV, V, new_add;
	for(int i = 0; i < nxt_gate_id.size(); ++i)
	{
		prime_field::field_element zero_value, one_value;
		int g = nxt_gate_id[i];
		int g_zero = g << 1, g_one = g << 1 | 1;
		if(mult_array.find(g_zero) != mult_array.end())
		{
			zero_value = mult_array[g_zero].eval(previous_random);
		}
		else
		{
			zero_value = prime_field::field_element(0);
		}
		if(mult_array.find(g_one) != mult_array.end())
		{
			one_value = mult_array[g_one].eval(previous_random);
		}
		else
		{
			one_value = prime_field::field_element(0);
		}
		mult_addV[g] = linear_poly(one_value - zero_value, zero_value);

		if(V_mult.find(g_zero) != V_mult.end())
		{
			zero_value = V_mult[g_zero].eval(previous_random);
		}
		else
		{
			zero_value = prime_field::field_element(0);
		}

		if(V_mult.find(g_one) != V_mult.end())
		{
			one_value = V_mult[g_one].eval(previous_random);
		}
		else
		{
			one_value = prime_field::field_element(0);
		}

		V[g] = linear_poly(one_value - zero_value, zero_value);
		ret = ret + mult_addV[g] * V[g];
	}
	//add gate

	mult_array = mult_addV;
	V_mult = V;
	mult_addV.clear(), V.clear();
	new_add.clear();
	for(int i = 0; i < nxt_gate_id.size(); ++i)
	{
		prime_field::field_element zero_value, one_value;
		int g = nxt_gate_id[i];
		int g_zero = g << 1, g_one = g << 1 | 1;
		if(addV_array.find(g_zero) != addV_array.end())
		{
			zero_value = addV_array[g_zero].eval(previous_random);
		}
		else
		{
			zero_value = prime_field::field_element(0);
		}
		if(addV_array.find(g_one) != addV_array.end())
		{
			one_value = addV_array[g_one].eval(previous_random);
		}
		else
		{
			one_value = prime_field::field_element(0);
		}
		mult_addV[g] = linear_poly(one_value - zero_value, zero_value);
		if(V_add.find(g_zero) != V_add.end())
		{
			zero_value = V_add[g_zero].eval(previous_random);
		}
		else
		{
			zero_value = prime_field::field_element(0);
		}
		if(V_add.find(g_one) != V_add.end())
		{
			one_value = V_add[g_one].eval(previous_random);
		}
		else
		{
			one_value = prime_field::field_element(0);
		}
		V[g] = linear_poly(one_value - zero_value, zero_value);

		if(add_array.find(g_zero) != add_array.end())
		{
			zero_value = add_array[g_zero].eval(previous_random);
		}
		else
		{
			zero_value = prime_field::field_element(0);
		}
		if(add_array.find(g_one) != add_array.end())
		{
			one_value = add_array[g_one].eval(previous_random);
		}
		else
		{
			one_value = prime_field::field_element(0);
		}
		new_add[g] = linear_poly(one_value - zero_value, zero_value);
		ret = ret + new_add[g] * V[g] + quadratic_poly(0, mult_addV[g].a, mult_addV[g].b);
	}
	addV_array = mult_addV;
	add_array = new_add;
	V_add = V;
	current_sumcheck_gates = nxt_gate_id;
	total_time += clock() - t0;
	return ret;
}
std::pair<prime_field::field_element, prime_field::field_element> prover::sumcheck_finalize(prime_field::field_element previous_random)
{
	v_v = (*V_add.begin()).second.eval(previous_random);
	return std::make_pair(v_u, v_v);
}
void prover::proof_init()
{
	//todo
}