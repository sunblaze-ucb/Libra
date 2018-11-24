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

void prover::init_array(int max_bit_length)
{
	mult_array = new linear_poly[(1 << max_bit_length)];
	add_array = new linear_poly[(1 << max_bit_length)];
	V_mult_add = new linear_poly[(1 << max_bit_length)];
	addV_array = new linear_poly[(1 << max_bit_length)];
	beta_g_r0 = new prime_field::field_element[(1 << max_bit_length)];
	beta_g_r1 = new prime_field::field_element[(1 << max_bit_length)];
	beta_g_sum = new prime_field::field_element[(1 << max_bit_length)];
	beta_u = new prime_field::field_element[(1 << max_bit_length)];
}

void prover::delete_self()
{
	delete[] mult_array;
	delete[] add_array;
	delete[] V_mult_add;
	delete[] addV_array;
	delete[] beta_g_r0;
	delete[] beta_g_r1;
	delete[] beta_u;
	delete[] beta_g_sum;
	for(int i = 0; i < C.total_depth; ++i)
		delete[] circuit_value[i];
}

prover::~prover()
{
}

void prover::sumcheck_phase1_init()
{
	fprintf(stderr, "sumcheck level %d, phase1 init start\n", sumcheck_layer_id);
	clock_t t0 = clock();
	//mult init
	clock_t t_init = clock();
	total_uv = (1 << C.circuit[sumcheck_layer_id - 1].bit_length);
	prime_field::field_element zero = prime_field::field_element(0);
	for(int i = 0; i < total_uv; ++i)
	{
		mult_array[i].a = zero;
		mult_array[i].b = zero;
		V_mult_add[i] = circuit_value[sumcheck_layer_id - 1][i];

		addV_array[i].a = zero;
		addV_array[i].b = zero;
		add_array[i].a = zero;
		add_array[i].b = zero;
	}
	fprintf(stderr, "	init time %f\n", ((float)clock() - t_init) / (float)CLOCKS_PER_SEC);
	clock_t t_beta = clock();
	beta_g_r0[0] = alpha;
	beta_g_r1[0] = beta;
	for(int i = 0; i < length_g; ++i)
	{
		for(int j = 0; j < (1 << i); ++j)
		{
			if(i % 4 == 0 || i + 1 == length_g)
			{
				beta_g_r0[j | (1 << i)].value = beta_g_r0[j].value * r_0[i].value % prime_field::mod;
				beta_g_r1[j | (1 << i)].value = beta_g_r1[j].value * r_1[i].value % prime_field::mod;
			}
			else
			{
				beta_g_r0[j | (1 << i)].value = beta_g_r0[j].value * r_0[i].value;
				beta_g_r1[j | (1 << i)].value = beta_g_r1[j].value * r_1[i].value;
			}
		}
		for(int j = 0; j < (1 << i); ++j)
		{
			if(i % 4 == 0 || i + 1 == length_g)
			{
				beta_g_r0[j].value = beta_g_r0[j].value * one_minus_r_0[i].value % prime_field::mod;
				beta_g_r1[j].value = beta_g_r1[j].value * one_minus_r_1[i].value % prime_field::mod;
			}
			else
			{
				beta_g_r0[j].value = beta_g_r0[j].value * one_minus_r_0[i].value;
				beta_g_r1[j].value = beta_g_r1[j].value * one_minus_r_1[i].value;
			}
		}
	}

	for(int i = 0; i < (1 << length_g); ++i)
	{
		beta_g_sum[i] = beta_g_r0[i] + beta_g_r1[i];
	}

	fprintf(stderr, "	beta calc time %f\n", ((float)clock() - t_beta) / (float)CLOCKS_PER_SEC);
	clock_t t_array = clock();
	for(int i = 0; i < (1 << length_g); ++i)
	{
		int u, v;
		u = C.circuit[sumcheck_layer_id].gates[i].u;
		v = C.circuit[sumcheck_layer_id].gates[i].v;
		if(C.circuit[sumcheck_layer_id].gates[i].ty == 0) //add gate
		{
			addV_array[u].b.value = addV_array[u].b.value + circuit_value[sumcheck_layer_id - 1][v].value * beta_g_sum[i].value;
			add_array[u].b.value = add_array[u].b.value + beta_g_sum[i].value;
		}
		if(C.circuit[sumcheck_layer_id].gates[i].ty == 1) //mult gate
		{
			mult_array[u].b.value = mult_array[u].b.value + circuit_value[sumcheck_layer_id - 1][v].value * beta_g_sum[i].value;
		}
	}
	fprintf(stderr, "	array calc time %f\n", ((float)clock() - t_array) / (float)CLOCKS_PER_SEC);

	fprintf(stderr, "sumcheck level %d, phase1 init finished\n", sumcheck_layer_id);
	clock_t time_e = clock() - t0;
	fprintf(stderr, "time %f\n", (float)time_e / (float)CLOCKS_PER_SEC);
	total_time += time_e;
}

quadratic_poly prover::sumcheck_phase1_update(prime_field::field_element previous_random, int current_bit)
{
	fprintf(stderr, "sumcheck level %d, phase1 update start\n", sumcheck_layer_id);
	clock_t t0 = clock();
	quadratic_poly ret = quadratic_poly(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));
	for(int i = 0; i < (total_uv >> 1); ++i)
	{
		prime_field::field_element zero_value, one_value;
		int g_zero = i << 1, g_one = i << 1 | 1;
		if(current_bit == 0)
		{
			mult_array[i].b = mult_array[g_zero].b;
			mult_array[i].a = mult_array[g_one].b - mult_array[i].b;
			
			V_mult_add[i].b = V_mult_add[g_zero].b;
			V_mult_add[i].a = V_mult_add[g_one].b - V_mult_add[i].b;

			addV_array[i].b = addV_array[g_zero].b;
			addV_array[i].a = addV_array[g_one].b - addV_array[i].b;
			add_array[i].b = add_array[g_zero].b;
			add_array[i].a = add_array[g_one].b - add_array[i].b;
		}
		else
		{
			mult_array[i].b = mult_array[g_zero].a * previous_random + mult_array[g_zero].b;
			mult_array[i].a = mult_array[g_one].a * previous_random + mult_array[g_one].b - mult_array[i].b;
			
			V_mult_add[i].b = V_mult_add[g_zero].a * previous_random + V_mult_add[g_zero].b;
			V_mult_add[i].a = V_mult_add[g_one].a * previous_random + V_mult_add[g_one].b - V_mult_add[i].b;

			addV_array[i].b = addV_array[g_zero].a * previous_random + addV_array[g_zero].b;
			addV_array[i].a = addV_array[g_one].a * previous_random + addV_array[g_one].b - addV_array[i].b;
			add_array[i].b = add_array[g_zero].a * previous_random + add_array[g_zero].b;
			add_array[i].a = add_array[g_one].a * previous_random + add_array[g_one].b - add_array[i].b;
		}
		ret.a.value = ret.a.value + (mult_array[i].a.value + add_array[i].a.value) * V_mult_add[i].a.value;
		ret.b.value = ret.b.value + (mult_array[i].a.value + add_array[i].a.value) * V_mult_add[i].b.value + (mult_array[i].b.value + add_array[i].b.value) * V_mult_add[i].a.value
								  + addV_array[i].a.value;
		ret.c.value = ret.c.value + (mult_array[i].b.value + add_array[i].b.value) * V_mult_add[i].b.value
								  + addV_array[i].b.value;
	}

	total_uv >>= 1;
	ret.a.value = ret.a.value % prime_field::mod;
	ret.b.value = ret.b.value % prime_field::mod;
	ret.c.value = ret.c.value % prime_field::mod;

	if(ret.a.value < prime_field::field_element(0).value)
		ret.a.value = ret.a.value + prime_field::mod;
	if(ret.b.value < prime_field::field_element(0).value)
		ret.b.value = ret.b.value + prime_field::mod;
	if(ret.c.value < prime_field::field_element(0).value)
		ret.c.value = ret.c.value + prime_field::mod;
	clock_t time_e = clock() - t0;
	fprintf(stderr, "sumcheck level %d, phase1 update finished\n", sumcheck_layer_id);
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

	//mult
	//init betag

	beta_u[0] = prime_field::field_element(1);
	for(int i = 0; i < length_u; ++i)
	{
		for(int j = 0; j < (1 << i); ++j)
			beta_u[j | (1 << i)] = beta_u[j] * r_u[i];
			
		for(int j = 0; j < (1 << i); ++j)
			beta_u[j] = beta_u[j] * one_minus_r_u[i];
	}
	
	total_uv = (1 << C.circuit[sumcheck_layer_id - 1].bit_length);
	int total_g = (1 << C.circuit[sumcheck_layer_id].bit_length);
	prime_field::field_element zero = prime_field::field_element(0);
	for(int i = 0; i < total_uv; ++i)
	{
		mult_array[i].a = zero;
		mult_array[i].b = zero;
		add_array[i].a = zero;
		add_array[i].b = zero;
		addV_array[i].a = zero;
		addV_array[i].b = zero;
		V_mult_add[i] = circuit_value[sumcheck_layer_id - 1][i];
	}

	for(int i = 0; i < total_g; ++i)
	{
		int ty = C.circuit[sumcheck_layer_id].gates[i].ty;
		int u = C.circuit[sumcheck_layer_id].gates[i].u;
		int v = C.circuit[sumcheck_layer_id].gates[i].v;
		if(ty == 1) //mult gate
		{
			mult_array[v].b.value = mult_array[v].b.value + (beta_g_r0[i].value + beta_g_r1[i].value) * beta_u[u].value * v_u.value;
		}
		if(ty == 0) //add gate
		{
			add_array[v].b.value = add_array[v].b.value + (beta_g_r0[i].value + beta_g_r1[i].value) * beta_u[u].value;
			addV_array[v].b.value = add_array[v].b.value * v_u.value + addV_array[v].b.value;
		}
	}

//	for(int i = 0; i < total_uv; ++i)
//	{
//		fprintf(stderr, "add %s\n", add_array[i].b.to_string(10).c_str());
//		fprintf(stderr, "addV %s\n", add_array[i].b.to_string(10).c_str());
//		fprintf(stderr, "V_mult_add %s\n", add_array[i].b.to_string(10).c_str());
//	}
	clock_t time_e = clock() - t0;
	fprintf(stderr, "sumcheck level %d, phase2 init finished\n", sumcheck_layer_id);
	fprintf(stderr, "time %f\n", (float)time_e / (float)CLOCKS_PER_SEC);
	total_time += time_e;
}

quadratic_poly prover::sumcheck_phase2_update(prime_field::field_element previous_random, int current_bit)
{
	fprintf(stderr, "sumcheck level %d, phase2 update start\n", sumcheck_layer_id);
	clock_t t0 = clock();
	quadratic_poly ret = quadratic_poly(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));
	for(int i = 0; i < (total_uv >> 1); ++i)
	{
		int g_zero = i << 1, g_one = i << 1 | 1;
		if(current_bit == 0)
		{
			mult_array[i].b = mult_array[g_zero].b;
			mult_array[i].a = mult_array[g_one].b - mult_array[i].b;
			
			V_mult_add[i].b = V_mult_add[g_zero].b;
			V_mult_add[i].a = V_mult_add[g_one].b - V_mult_add[i].b;

			addV_array[i].b = addV_array[g_zero].b;
			addV_array[i].a = addV_array[g_one].b - addV_array[i].b;

			add_array[i].b = add_array[g_zero].b;
			add_array[i].a = add_array[g_one].b - add_array[i].b;
		}
		else
		{
			mult_array[i].b = mult_array[g_zero].a * previous_random + mult_array[g_zero].b;
			mult_array[i].a = mult_array[g_one].a * previous_random + mult_array[g_one].b - mult_array[i].b;
			
			V_mult_add[i].b = V_mult_add[g_zero].a * previous_random + V_mult_add[g_zero].b;
			V_mult_add[i].a = V_mult_add[g_one].a * previous_random + V_mult_add[g_one].b - V_mult_add[i].b;

			addV_array[i].b = addV_array[g_zero].a * previous_random + addV_array[g_zero].b;
			addV_array[i].a = addV_array[g_one].a * previous_random + addV_array[g_one].b - addV_array[i].b;

			add_array[i].b = add_array[g_zero].a * previous_random + add_array[g_zero].b;
			add_array[i].a = add_array[g_one].a * previous_random + add_array[g_one].b - add_array[i].b;
		}
		ret.a.value = ret.a.value + (mult_array[i].a.value + add_array[i].a.value) * V_mult_add[i].a.value;
		ret.b.value = ret.b.value + (mult_array[i].a.value + add_array[i].a.value) * V_mult_add[i].b.value
								  +	(mult_array[i].b.value + add_array[i].b.value) * V_mult_add[i].a.value
								  + addV_array[i].a.value;
		ret.c.value = ret.c.value + (mult_array[i].b.value + add_array[i].b.value) * V_mult_add[i].b.value
								  + addV_array[i].b.value;
	}

	total_uv >>= 1;
	ret.a.value = ret.a.value % prime_field::mod;
	ret.b.value = ret.b.value % prime_field::mod;
	ret.c.value = ret.c.value % prime_field::mod;
	if(ret.a.value < prime_field::field_element(0).value)
		ret.a.value = ret.a.value + prime_field::mod;
	if(ret.b.value < prime_field::field_element(0).value)
		ret.b.value = ret.b.value + prime_field::mod;
	if(ret.c.value < prime_field::field_element(0).value)
		ret.c.value = ret.c.value + prime_field::mod;
	clock_t time_e = clock() - t0;
	fprintf(stderr, "sumcheck level %d, phase2 update finished\n", sumcheck_layer_id);
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
