#include "linear_gkr/prover_fast_track.h"

void prover::get_circuit(const layered_circuit &from_verifier)
{
	C = from_verifier;
}

prime_field::field_element prover::V_res(const prime_field::field_element* one_minus_r_0, const prime_field::field_element* r_0, const prime_field::field_element* output_raw, int r_0_size, int output_size)
{
	clock_t t0 = clock();
	prime_field::field_element *output;
	output = new prime_field::field_element[output_size];
	for(int i = 0; i < output_size; ++i)
		output[i] = output_raw[i];
	for(int i = 0; i < r_0_size; ++i)
	{
		for(int j = 0; j < (output_size >> 1); ++j)
		{
			output[j] = output[j << 1] * one_minus_r_0[i] + output[j << 1 | 1] * r_0[i];
		}
	}
	total_time += (clock() - t0);
	prime_field::field_element res = output[0];
	delete[] output;
	return res;
}

prime_field::field_element* prover::evaluate()
{
	clock_t t0 = clock();
	circuit_value[0] = new prime_field::field_element[(1 << C.circuit[0].bit_length)];
	for(int i = 0; i < (1 << C.circuit[0].bit_length); ++i)
	{
		int g, u, v, ty;
		g = i;
		u = C.circuit[0].gates[g].u;
		v = C.circuit[0].gates[g].v;
		ty = C.circuit[0].gates[g].ty;
		assert(ty == 3);
		circuit_value[0][g] = mpz_class(u);
	}
	for(int i = 1; i < C.total_depth; ++i)
	{
		circuit_value[i] = new prime_field::field_element[(1 << C.circuit[i].bit_length)];
		for(int j = 0; j < (1 << C.circuit[i].bit_length); ++j)
		{
			int g, u, v, ty;
			g = j;
			ty = C.circuit[i].gates[g].ty;
			u = C.circuit[i].gates[g].u;
			v = C.circuit[i].gates[g].v;
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
				circuit_value[i][g] = prime_field::field_element(0);
			}
			else if(ty == 3)
			{
				circuit_value[i][g] = prime_field::field_element(u);
			}
			else
			{
				assert(false);
			}
		}
	}
/*
	//cheating prover test
	ret[0].second =  prime_field::field_element(111);
*/
	fprintf(stderr, "total evaluation time: %f\n", ((float)clock() - (float)t0) / CLOCKS_PER_SEC);
	return circuit_value[C.total_depth - 1];
}

void prover::sumcheck_init(int layer_id, int bit_length_g, int bit_length_u, int bit_length_v, 
	const prime_field::field_element &a, const prime_field::field_element &b, 
	const prime_field::field_element* R_0, const prime_field::field_element* R_1,
	prime_field::field_element* o_r_0, prime_field::field_element *o_r_1)
{
	r_0 = R_0;
	r_1 = R_1;
	alpha = a;
	beta = b;
	sumcheck_layer_id = layer_id;
	length_g = bit_length_g;
	length_u = bit_length_u;
	length_v = bit_length_v;
	one_minus_r_0 = o_r_0;
	one_minus_r_1 = o_r_1;
}
/*
void prover::DFS(linear_poly* AddV, linear_poly *Mult, linear_poly *Add, int g, int depth, 
	prime_field::field_element alpha_value, prime_field::field_element beta_value)
{
	if(depth == length_g)
	{
		int u, v;
		u = C.circuit[sumcheck_layer_id].gates[g].u;
		v = C.circuit[sumcheck_layer_id].gates[g].v;
		if(C.circuit[sumcheck_layer_id].gates[g].ty == 0) //add gate
		{
			AddV[u] = AddV[u] + circuit_value[sumcheck_layer_id - 1][v] * (alpha_value * alpha + beta_value * beta);
			Add[u] = Add[u] + (alpha_value * alpha + beta_value * beta);
		}
		else if(C.circuit[sumcheck_layer_id].gates[g].ty == 1) //mult gate
		{
			Mult[u] = Mult[u] + circuit_value[sumcheck_layer_id - 1][v] * (alpha_value * alpha + beta_value * beta);
		}
	}
	else
	{
		g &= ((-1) ^ (1 << depth));
		DFS(AddV, Mult, Add, g, depth + 1, 
			alpha_value * (one_minus_r_0[depth]), beta_value * (one_minus_r_1[depth]));
		g |= (1 << depth);
		DFS(AddV, Mult, Add, g, depth + 1, alpha_value * r_0[depth], beta_value * r_1[depth]);
	}
}
*/

void prover::init_array(int max_bit_length)
{
	mult_array = new linear_poly[(1 << max_bit_length)];
	add_array = new linear_poly[(1 << max_bit_length)];
	V_mult_add = new linear_poly[(1 << max_bit_length)];
	addV_array = new linear_poly[(1 << max_bit_length)];
	beta_g_r0 = new prime_field::field_element[(1 << max_bit_length)];
	beta_g_r1 = new prime_field::field_element[(1 << max_bit_length)];
	beta_u = new prime_field::field_element[(1 << max_bit_length)];
}

prover::~prover()
{
	delete[] mult_array;
	delete[] add_array;
	delete[] V_mult_add;
	delete[] addV_array;
	delete[] beta_g_r0;
	delete[] beta_g_r1;
	delete[] beta_u;
	for(int i = 0; i < C.total_depth; ++i)
		delete[] circuit_value[i];
}

void prover::sumcheck_phase1_init()
{
	fprintf(stderr, "sumcheck level %d, phase1 init start\n", sumcheck_layer_id);
	clock_t t0 = clock();
	//mult init
	total_uv = (1 << C.circuit[sumcheck_layer_id - 1].bit_length);
	for(int i = 0; i < total_uv; ++i)
	{
		int g = i;
		mult_array[g] = 
			linear_poly(prime_field::field_element(0), prime_field::field_element(0));
		V_mult_add[g] = 
			circuit_value[sumcheck_layer_id - 1][g];
	}
	//add init

	for(int i = 0; i < total_uv; ++i)
	{
		int g = i;
		addV_array[g] = 
			linear_poly(prime_field::field_element(0), prime_field::field_element(0));
		add_array[g] = 
			linear_poly(prime_field::field_element(0), prime_field::field_element(0));
	}

	//DFS(addV_array, mult_array, add_array, 0, 0, prime_field::field_element(1), prime_field::field_element(1));
	
	beta_g_r0[0] = prime_field::field_element(1);
	beta_g_r1[0] = prime_field::field_element(1);
	for(int i = 0; i < length_g; ++i)
	{
		for(int j = 0; j < (1 << i); ++j)
		{
			beta_g_r0[j | (1 << i)] = beta_g_r0[j] * r_0[i];
			beta_g_r0[j] = beta_g_r0[j] * one_minus_r_0[i];
			beta_g_r1[j | (1 << i)] = beta_g_r1[j] * r_1[i];
			beta_g_r1[j] = beta_g_r1[j] * one_minus_r_1[i];
		}
	}

	for(int i = 0; i < (1 << length_g); ++i)
	{
		int u, v;
		u = C.circuit[sumcheck_layer_id].gates[i].u;
		v = C.circuit[sumcheck_layer_id].gates[i].v;
		if(C.circuit[sumcheck_layer_id].gates[i].ty == 0) //add gate
		{
			addV_array[u] = addV_array[u] + circuit_value[sumcheck_layer_id - 1][v] * (beta_g_r0[i] * alpha + beta_g_r1[i] * beta);
			add_array[u] = add_array[u] + (beta_g_r0[i] * alpha + beta_g_r1[i] * beta);
		}
		else if(C.circuit[sumcheck_layer_id].gates[i].ty == 1) //mult gate
		{
			mult_array[u] = mult_array[u] + circuit_value[sumcheck_layer_id - 1][v] * (beta_g_r0[i] * alpha + beta_g_r1[i] * beta);
		}
	}

	fprintf(stderr, "sumcheck level %d, phase1 init finished\n", sumcheck_layer_id);
	clock_t time_e = clock() - t0;
	fprintf(stderr, "time %f\n", (float)time_e / (float)CLOCKS_PER_SEC);
	total_time += time_e;
}

quadratic_poly prover::sumcheck_phase1_update(prime_field::field_element previous_random)
{
	clock_t t0 = clock();
	quadratic_poly ret = quadratic_poly(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));
	for(int i = 0; i < (total_uv >> 1); ++i)
	{
		prime_field::field_element zero_value, one_value;
		int g = i;
		int g_zero = g << 1, g_one = g << 1 | 1;
		zero_value = mult_array[g_zero].a * previous_random + mult_array[g_zero].b;
		one_value = mult_array[g_one].a * previous_random + mult_array[g_one].b;
		mult_array[g] = linear_poly(one_value - zero_value, zero_value);

		zero_value = V_mult_add[g_zero].a * previous_random + V_mult_add[g_zero].b;
		one_value = V_mult_add[g_one].a * previous_random + V_mult_add[g_one].b;
		V_mult_add[g] = linear_poly(one_value - zero_value, zero_value);
		//ret = ret + mult_array[g] * V_mult_add[g];

		ret.a.value = ret.a.value + mult_array[g].a.value * V_mult_add[g].a.value;
		ret.b.value = ret.b.value + mult_array[g].a.value * V_mult_add[g].b.value + mult_array[g].b.value * V_mult_add[g].a.value;
		ret.c.value = ret.c.value + mult_array[g].b.value * V_mult_add[g].b.value;
	}

	//add gate
	for(int i = 0; i < (total_uv >> 1); ++i)
	{
		prime_field::field_element zero_value, one_value;
		int g = i;
		int g_zero = g << 1, g_one = g << 1 | 1;
		zero_value = addV_array[g_zero].a * previous_random + addV_array[g_zero].b;
		one_value = addV_array[g_one].a * previous_random + addV_array[g_one].b;
		addV_array[g] = linear_poly(one_value - zero_value, zero_value);
		zero_value = add_array[g_zero].a * previous_random + add_array[g_zero].b;
		one_value = add_array[g_one].a * previous_random + add_array[g_one].b;
		add_array[g] = linear_poly(one_value - zero_value, zero_value);
		//ret = ret + add_array[g] * V_mult_add[g] + quadratic_poly(0, addV_array[g].a, addV_array[g].b);
		ret.a.value = ret.a.value + add_array[g].a.value * V_mult_add[g].a.value;
		ret.b.value = ret.b.value + add_array[g].a.value * V_mult_add[g].b.value + add_array[g].b.value * V_mult_add[g].a.value;
		ret.c.value = ret.c.value + add_array[g].b.value * V_mult_add[g].b.value;

		ret.c.value = ret.c.value + addV_array[g].b.value;
		ret.b.value = ret.b.value + addV_array[g].a.value;
	}
	total_uv >>= 1;
	ret.a.value = ret.a.value % prime_field::mod;
	ret.b.value = ret.b.value % prime_field::mod;
	ret.c.value = ret.c.value % prime_field::mod;
	clock_t time_e = clock() - t0;
	fprintf(stderr, "time %f\n", (float)time_e / (float)CLOCKS_PER_SEC);
	total_time += time_e;
	return ret;
}

void prover::sumcheck_phase2_init(prime_field::field_element previous_random, const prime_field::field_element* r_u, const prime_field::field_element* one_minus_r_u)
{
	fprintf(stderr, "sumcheck level %d, phase2 init start\n", sumcheck_layer_id);
	clock_t t0 = clock();
	v_u = V_mult_add[0].eval(previous_random);
//	fprintf(stderr, "v_u %s\n", v_u.to_string(10).c_str());
	
//	for(int i = 0; i < length_g; ++i)
//	{
//		if(r_0[i] != (prime_field::field_element(1) - one_minus_r_0[i]))
//		{
//			assert(false);
//		}
//		if(r_1[i] != (prime_field::field_element(1) - one_minus_r_1[i]))
//		{
//			assert(false);
//		}
//	}
	//mult
	//init betag

//	beta_g_r0[0] = prime_field::field_element(1);
//	beta_g_r1[0] = prime_field::field_element(1);

//	for(int i = 0; i < length_g; ++i)
//	{
//		for(int j = 0; j < (1 << i); ++j)
//		{
//			beta_g_r0[j | (1 << i)] = beta_g_r0[j] * r_0[i];
//			beta_g_r0[j] = beta_g_r0[j] * one_minus_r_0[i];
//			beta_g_r1[j | (1 << i)] = beta_g_r1[j] * r_1[i];
//			beta_g_r1[j] = beta_g_r1[j] * one_minus_r_1[i];
//		}
//	}
	beta_u[0] = prime_field::field_element(1);
	for(int i = 0; i < length_u; ++i)
	{
		for(int j = 0; j < (1 << i); ++j)
		{
			beta_u[j | (1 << i)] = beta_u[j] * r_u[i];
			beta_u[j] = beta_u[j] * one_minus_r_u[i];
		}
	}
	
	//DFS_betag(beta_g_r0, r_0, one_minus_r_0, 0, 0, prime_field::field_element(1), 1);
	//DFS_betag(beta_g_r1, r_1, one_minus_r_1, 0, 0, prime_field::field_element(1), 1);
	//init betau
	//DFS_betau(beta_u, r_u, one_minus_r_u, 0, 0, prime_field::field_element(1), 1);

	total_uv = (1 << C.circuit[sumcheck_layer_id - 1].bit_length);
	int total_g = (1 << C.circuit[sumcheck_layer_id].bit_length);

	for(int i = 0; i < total_uv; ++i)
	{
		int v = i;
		mult_array[v] = linear_poly(prime_field::field_element(0), prime_field::field_element(0));
		add_array[v] = linear_poly(prime_field::field_element(0), prime_field::field_element(0));
		addV_array[v] = linear_poly(prime_field::field_element(0), prime_field::field_element(0));
		V_mult_add[v] = circuit_value[sumcheck_layer_id - 1][v];
	}

	for(int i = 0; i < total_g; ++i)
	{
		int g = i;
		int ty = C.circuit[sumcheck_layer_id].gates[g].ty;
		int u = C.circuit[sumcheck_layer_id].gates[g].u;
		int v = C.circuit[sumcheck_layer_id].gates[g].v;
		if(ty == 1) //mult gate
		{
			//mult_array[v] = mult_array[v] + linear_poly((beta_g_r0[g] * beta_u[u] * alpha + beta_g_r1[g] * beta_u[u] * beta) * v_u);
			mult_array[v].b.value = mult_array[v].b.value + (beta_g_r0[g].value * beta_u[u].value * alpha.value + beta_g_r1[g].value * beta_u[u].value * beta.value) * v_u.value;
		}
		if(ty == 0) //add gate
		{
		//	add_array[v] = add_array[v] + linear_poly(beta_g_r0[g] * beta_u[u] * alpha + beta_g_r1[g] * beta_u[u] * beta);
		//	addV_array[v] = addV_array[v] + linear_poly((beta_g_r0[g] * beta_u[u] * alpha + beta_g_r1[g] * beta_u[u] * beta) * v_u);
			add_array[v].b.value = add_array[v].b.value + beta_g_r0[g].value * beta_u[u].value * alpha.value + beta_g_r1[g].value * beta_u[u].value * beta.value;
			addV_array[v].b.value = add_array[v].b.value * v_u.value + addV_array[v].b.value;
		}
	}
//	DFS_betag(beta_g_r0, r_0, one_minus_r_0, 0, 0, prime_field::field_element(1), 0);
//	DFS_betag(beta_g_r1, r_1, one_minus_r_1, 0, 0, prime_field::field_element(1), 0);
//	DFS_betau(beta_u, r_u, one_minus_r_u, 0, 0, prime_field::field_element(1), 0);

	for(int i = 0; i < total_uv; ++i)
	{
		int v = i;
		mult_array[v].b.value = mult_array[v].b.value % prime_field::mod;
		add_array[v].b.value = add_array[v].b.value % prime_field::mod;
		addV_array[v].b.value = addV_array[v].b.value % prime_field::mod;
	}

//	for(int i = 0; i < total_uv; ++i)
//	{
//		fprintf(stderr, "add %s\n", add_array[i].b.to_string(10).c_str());
//		fprintf(stderr, "addV %s\n", add_array[i].b.to_string(10).c_str());
//		fprintf(stderr, "V_mult_add %s\n", add_array[i].b.to_string(10).c_str());
//	}
	clock_t time_e = clock() - t0;
	fprintf(stderr, "time %f\n", (float)time_e / (float)CLOCKS_PER_SEC);
	total_time += time_e;
}

quadratic_poly prover::sumcheck_phase2_update(prime_field::field_element previous_random)
{
	clock_t t0 = clock();
	quadratic_poly ret = quadratic_poly(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));
	for(int i = 0; i < (total_uv >> 1); ++i)
	{
		prime_field::field_element zero_value, one_value;
		int g = i;
		int g_zero = g << 1, g_one = g << 1 | 1;
		zero_value = mult_array[g_zero].a * previous_random + mult_array[g_zero].b;
		one_value = mult_array[g_one].a * previous_random + mult_array[g_one].b;
		mult_array[g] = linear_poly(one_value - zero_value, zero_value);
		
		zero_value = V_mult_add[g_zero].a * previous_random + V_mult_add[g_zero].b;
		one_value = V_mult_add[g_one].a * previous_random + V_mult_add[g_one].b;
		V_mult_add[g] = linear_poly(one_value - zero_value, zero_value);
		
		ret.a.value = ret.a.value + mult_array[g].a.value * V_mult_add[g].a.value;
		ret.b.value = ret.b.value + mult_array[g].a.value * V_mult_add[g].b.value + 
									mult_array[g].b.value * V_mult_add[g].a.value;
		ret.c.value = ret.c.value + mult_array[g].b.value * V_mult_add[g].b.value;
		//ret = ret + mult_array[g] * V_mult_add[g];
	}
	//add gate

	for(int i = 0; i < (total_uv >> 1); ++i)
	{
		prime_field::field_element zero_value, one_value;
		int g = i;
		int g_zero = g << 1, g_one = g << 1 | 1;
		zero_value = addV_array[g_zero].a * previous_random + addV_array[g_zero].b;
		one_value = addV_array[g_one].a * previous_random + addV_array[g_one].b;
		addV_array[g] = linear_poly(one_value - zero_value, zero_value);

		zero_value = add_array[g_zero].a * previous_random + add_array[g_zero].b;
		one_value = add_array[g_one].a * previous_random + add_array[g_one].b;
		add_array[g] = linear_poly(one_value - zero_value, zero_value);

		ret.a.value = ret.a.value + add_array[g].a.value * V_mult_add[g].a.value;
		ret.b.value = ret.b.value + add_array[g].a.value * V_mult_add[g].b.value + 
									add_array[g].b.value * V_mult_add[g].a.value;
		ret.c.value = ret.c.value + add_array[g].b.value * V_mult_add[g].b.value;
		ret.b.value = ret.b.value + addV_array[g].a.value;
		ret.c.value = ret.c.value + addV_array[g].b.value;
		//ret = ret + add_array[g] * V_mult_add[g] + quadratic_poly(0, addV_array[g].a, addV_array[g].b);
	}
	total_uv >>= 1;
	ret.a.value = ret.a.value % prime_field::mod;
	ret.b.value = ret.b.value % prime_field::mod;
	ret.c.value = ret.c.value % prime_field::mod;
	clock_t time_e = clock() - t0;
	fprintf(stderr, "time %f\n", (float)time_e / (float)CLOCKS_PER_SEC);
	total_time += time_e;
	return ret;
}
std::pair<prime_field::field_element, prime_field::field_element> prover::sumcheck_finalize(prime_field::field_element previous_random)
{
	v_v = V_mult_add[0].eval(previous_random);
	return std::make_pair(v_u, v_v);
}
void prover::proof_init()
{
	//todo
}