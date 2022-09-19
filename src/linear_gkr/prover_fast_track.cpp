#include "linear_gkr/prover_fast_track.h"

void prover::get_circuit(const layered_circuit &from_verifier)
{
	C = from_verifier;
}

prime_field::field_element prover::V_res(const prime_field::field_element* one_minus_r_0, const prime_field::field_element* r_0, const prime_field::field_element* output_raw, int r_0_size, int output_size)
{
	std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
	prime_field::field_element *output;
	output = new prime_field::field_element[output_size];
	for(int i = 0; i < output_size; ++i)
		output[i] = output_raw[i];
	for(int i = 0; i < r_0_size; ++i)
	{
		for(int j = 0; j < (output_size >> 1); ++j)
		{
			output[j].value = (output[j << 1].value * one_minus_r_0[i].value + output[j << 1 | 1].value * r_0[i].value) % prime_field::mod;
		}
		output_size >>= 1;
	}
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	total_time += time_span.count();
	prime_field::field_element res = output[0];
	delete[] output;

	if(res.value < 0)
		res.value = res.value + prime_field::mod;
	return res;
}

prime_field::field_element* prover::evaluate()
{
	std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
	circuit_value[0] = new prime_field::field_element[(1 << C.circuit[0].bit_length)];
	for(int i = 0; i < (1 << C.circuit[0].bit_length); ++i)
	{
		int g, u, v, ty;
		g = i;
		u = C.circuit[0].gates[g].u;
		v = C.circuit[0].gates[g].v;
		ty = C.circuit[0].gates[g].ty;
		assert(ty == 3 || ty == 2);
		circuit_value[0][g] = prime_field::field_element(u);
	}
	assert(C.total_depth < 1000000);
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
			else if(ty == 4)
			{
				circuit_value[i][g] = circuit_value[i - 1][u];
			}
			else if(ty == 5)
			{
				circuit_value[i][g] = prime_field::field_element(0);
				for(int k = u; k < v; ++k)
					circuit_value[i][g] = circuit_value[i][g] + circuit_value[i - 1][k];
			}
			else if(ty == 6)
			{
				circuit_value[i][g] = prime_field::field_element(1) - circuit_value[i - 1][u];
			}
			else if(ty == 7)
			{
				circuit_value[i][g] = circuit_value[i - 1][u] - circuit_value[i - 1][v];
			}
			else if(ty == 8)
			{
				auto &x = circuit_value[i - 1][u], &y = circuit_value[i - 1][v];
				circuit_value[i][g] = x + y - prime_field::field_element(2) * x * y;
			}
			else if(ty == 9)
			{
				auto &x = circuit_value[i - 1][u], &y = circuit_value[i - 1][v];
				circuit_value[i][g] = y - x * y;
			}
			else if(ty == 10)
			{
				circuit_value[i][g] = circuit_value[i - 1][u];
			}
			else
			{
				assert(false);
			}
		}
	}

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	std::cerr << "total evaluation time: " << time_span.count() << " seconds." << std::endl;
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
	add_mult_sum = new linear_poly[(1 << max_bit_length)];
	V_mult_add = new linear_poly[(1 << max_bit_length)];
	addV_array = new linear_poly[(1 << max_bit_length)];
	beta_g_sum = new prime_field::field_element[(1 << max_bit_length)];
	beta_u = new prime_field::field_element[(1 << max_bit_length)];
	int half_length = (max_bit_length >> 1) + 1;
	beta_g_r0_fhalf = new prime_field::field_element[(1 << half_length)];
	beta_g_r0_shalf = new prime_field::field_element[(1 << half_length)];
	beta_g_r1_fhalf = new prime_field::field_element[(1 << half_length)];
	beta_g_r1_shalf = new prime_field::field_element[(1 << half_length)];
	beta_u_fhalf = new prime_field::field_element[(1 << half_length)];
	beta_u_shalf = new prime_field::field_element[(1 << half_length)];
}

void prover::delete_self()
{
	delete[] add_mult_sum;
	delete[] V_mult_add;
	delete[] addV_array;
	delete[] beta_u;
	delete[] beta_g_sum;

	delete[] beta_g_r0_fhalf;
	delete[] beta_g_r0_shalf;
	delete[] beta_g_r1_fhalf;
	delete[] beta_g_r1_shalf;
	delete[] beta_u_fhalf;
	delete[] beta_u_shalf;
	for(int i = 0; i < C.total_depth; ++i)
		delete[] circuit_value[i];
}

prover::~prover()
{
}

void prover::sumcheck_phase1_init()
{
	std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
	//mult init
	total_uv = (1 << C.circuit[sumcheck_layer_id - 1].bit_length);
	prime_field::field_element zero = prime_field::field_element(0);
	for(int i = 0; i < total_uv; ++i)
	{
		V_mult_add[i] = circuit_value[sumcheck_layer_id - 1][i];

		addV_array[i].a = zero;
		addV_array[i].b = zero;
		add_mult_sum[i].a = zero;
		add_mult_sum[i].b = zero;
	}
	
	beta_g_r0_fhalf[0] = alpha;
	beta_g_r1_fhalf[0] = beta;
	beta_g_r0_shalf[0] = prime_field::field_element(1);
	beta_g_r1_shalf[0] = prime_field::field_element(1);

	int first_half = length_g >> 1, second_half = length_g - first_half;

	for(int i = 0; i < first_half; ++i)
	{
		for(int j = 0; j < (1 << i); ++j)
		{
			beta_g_r0_fhalf[j | (1 << i)].value = beta_g_r0_fhalf[j].value * r_0[i].value % prime_field::mod;
			beta_g_r0_fhalf[j].value = beta_g_r0_fhalf[j].value * one_minus_r_0[i].value % prime_field::mod;
			beta_g_r1_fhalf[j | (1 << i)].value = beta_g_r1_fhalf[j].value * r_1[i].value % prime_field::mod;
			beta_g_r1_fhalf[j].value = beta_g_r1_fhalf[j].value * one_minus_r_1[i].value % prime_field::mod;
		}
	}

	for(int i = 0; i < second_half; ++i)
	{
		for(int j = 0; j < (1 << i); ++j)
		{
			beta_g_r0_shalf[j | (1 << i)].value = beta_g_r0_shalf[j].value * r_0[i + first_half].value % prime_field::mod;
			beta_g_r0_shalf[j].value = beta_g_r0_shalf[j].value * one_minus_r_0[i + first_half].value % prime_field::mod;
			beta_g_r1_shalf[j | (1 << i)].value = beta_g_r1_shalf[j].value * r_1[i + first_half].value % prime_field::mod;
			beta_g_r1_shalf[j].value = beta_g_r1_shalf[j].value * one_minus_r_1[i + first_half].value % prime_field::mod;
		}
	}

	int mask_fhalf = (1 << first_half) - 1;

	for(int i = 0; i < (1 << length_g); ++i)
	{
		beta_g_sum[i].value = (beta_g_r0_fhalf[i & mask_fhalf].value * beta_g_r0_shalf[i >> first_half].value 
							 + beta_g_r1_fhalf[i & mask_fhalf].value * beta_g_r1_shalf[i >> first_half].value) % prime_field::mod;
	}
	for(int i = 0; i < (1 << length_g); ++i)
	{
		int u, v;
		u = C.circuit[sumcheck_layer_id].gates[i].u;
		v = C.circuit[sumcheck_layer_id].gates[i].v;
		switch(C.circuit[sumcheck_layer_id].gates[i].ty)
		{
			case 0: //add gate
			{
				auto tmp = (beta_g_r0_fhalf[i & mask_fhalf].value * beta_g_r0_shalf[i >> first_half].value 
						+ beta_g_r1_fhalf[i & mask_fhalf].value * beta_g_r1_shalf[i >> first_half].value) % prime_field::mod;
				addV_array[u].b.value = (addV_array[u].b.value + circuit_value[sumcheck_layer_id - 1][v].value * tmp) % prime_field::mod;
				add_mult_sum[u].b.value = (add_mult_sum[u].b.value + tmp) % prime_field::mod;
				break;
			}
			case 1: //mult gate
			{
				auto tmp = (beta_g_r0_fhalf[i & mask_fhalf].value * beta_g_r0_shalf[i >> first_half].value 
						+ beta_g_r1_fhalf[i & mask_fhalf].value * beta_g_r1_shalf[i >> first_half].value) % prime_field::mod;
				add_mult_sum[u].b.value = (add_mult_sum[u].b.value + circuit_value[sumcheck_layer_id - 1][v].value * tmp) % prime_field::mod;
				break;
			}
			case 4: //direct relay gate
			{
				auto tmp = (beta_g_r0_fhalf[u & mask_fhalf].value * beta_g_r0_shalf[u >> first_half].value 
						+ beta_g_r1_fhalf[u & mask_fhalf].value * beta_g_r1_shalf[u >> first_half].value) % prime_field::mod;
				add_mult_sum[u].b.value = (add_mult_sum[u].b.value + tmp) % prime_field::mod;
				break;
			}
			case 5: //sum gate
			{
				auto tmp = beta_g_r0_fhalf[i & mask_fhalf].value * beta_g_r0_shalf[i >> first_half].value 
					+ beta_g_r1_fhalf[i & mask_fhalf].value * beta_g_r1_shalf[i >> first_half].value;
				tmp = tmp % prime_field::mod;
				for(int j = u; j < v; ++j)
				{
					add_mult_sum[j].b.value = (add_mult_sum[j].b.value + tmp);
					if(add_mult_sum[j].b.value >= prime_field::mod)
						add_mult_sum[j].b.value = add_mult_sum[j].b.value - prime_field::mod;
				}
				break;
			}
			case 6: //NOT gate
			{
				auto tmp = (beta_g_r0_fhalf[i & mask_fhalf].value * beta_g_r0_shalf[i >> first_half].value 
						+ beta_g_r1_fhalf[i & mask_fhalf].value * beta_g_r1_shalf[i >> first_half].value) % prime_field::mod;
				add_mult_sum[u].b.value = (add_mult_sum[u].b.value + prime_field::mod - tmp) % prime_field::mod;
				addV_array[u].b.value = (addV_array[u].b.value + tmp) % prime_field::mod;
				break;
			}
			case 7: //minus gate
			{
				auto tmp = (beta_g_r0_fhalf[i & mask_fhalf].value * beta_g_r0_shalf[i >> first_half].value 
						+ beta_g_r1_fhalf[i & mask_fhalf].value * beta_g_r1_shalf[i >> first_half].value) % prime_field::mod;
				addV_array[u].b.value = (addV_array[u].b.value + prime_field::mod - (circuit_value[sumcheck_layer_id - 1][v].value * tmp % prime_field::mod)) % prime_field::mod;
				add_mult_sum[u].b.value = (add_mult_sum[u].b.value + tmp) % prime_field::mod;
				break;
			}
			case 8: //XOR gate
			{
				auto tmp = (beta_g_r0_fhalf[i & mask_fhalf].value * beta_g_r0_shalf[i >> first_half].value 
						+ beta_g_r1_fhalf[i & mask_fhalf].value * beta_g_r1_shalf[i >> first_half].value) % prime_field::mod;
				addV_array[u].b.value = (addV_array[u].b.value + tmp * circuit_value[sumcheck_layer_id - 1][v].value % prime_field::mod) % prime_field::mod;
				add_mult_sum[u].b.value = (add_mult_sum[u].b.value + tmp + prime_field::mod - (prime_field::field_element(2).value * circuit_value[sumcheck_layer_id - 1][v].value % prime_field::mod * tmp % prime_field::mod)) % prime_field::mod;
				break;
			}
			case 9: //NAAB gate
			{
				auto tmp = (beta_g_r0_fhalf[i & mask_fhalf].value * beta_g_r0_shalf[i >> first_half].value 
						+ beta_g_r1_fhalf[i & mask_fhalf].value * beta_g_r1_shalf[i >> first_half].value) % prime_field::mod;
				addV_array[u].b.value = (addV_array[u].b.value + tmp * circuit_value[sumcheck_layer_id - 1][v].value % prime_field::mod) % prime_field::mod;
				add_mult_sum[u].b.value = (add_mult_sum[u].b.value + prime_field::mod - (circuit_value[sumcheck_layer_id - 1][v].value * tmp % prime_field::mod)) % prime_field::mod;
				break;
			}
			case 10: //relay gate
			{
				auto tmp = (beta_g_r0_fhalf[i & mask_fhalf].value * beta_g_r0_shalf[i >> first_half].value 
						+ beta_g_r1_fhalf[i & mask_fhalf].value * beta_g_r1_shalf[i >> first_half].value) % prime_field::mod;
				add_mult_sum[u].b.value = (add_mult_sum[u].b.value + tmp) % prime_field::mod;
				break;
			}
			default:
			{
				break;
			}
		}
	}
	
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	double tmp = time_span.count();
	total_time += time_span.count();
}

quadratic_poly prover::sumcheck_phase1_update(prime_field::field_element previous_random, int current_bit)
{
	std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
	quadratic_poly ret = quadratic_poly(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));
	for(int i = 0; i < (total_uv >> 1); ++i)
	{
		prime_field::field_element zero_value, one_value;
		int g_zero = i << 1, g_one = i << 1 | 1;
		if(current_bit == 0)
		{
			V_mult_add[i].b = V_mult_add[g_zero].b;
			V_mult_add[i].a = V_mult_add[g_one].b - V_mult_add[i].b;

			addV_array[i].b = addV_array[g_zero].b;
			addV_array[i].a = addV_array[g_one].b - addV_array[i].b;

			add_mult_sum[i].b = add_mult_sum[g_zero].b;
			add_mult_sum[i].a = add_mult_sum[g_one].b - add_mult_sum[i].b;

		}
		else
		{
			V_mult_add[i].b.value = (V_mult_add[g_zero].a.value * previous_random.value + V_mult_add[g_zero].b.value) % prime_field::mod;
			V_mult_add[i].a.value = (V_mult_add[g_one].a.value * previous_random.value + V_mult_add[g_one].b.value - V_mult_add[i].b.value + prime_field::mod) % prime_field::mod;

			addV_array[i].b.value = (addV_array[g_zero].a.value * previous_random.value + addV_array[g_zero].b.value) % prime_field::mod;
			addV_array[i].a.value = (addV_array[g_one].a.value * previous_random.value + addV_array[g_one].b.value - addV_array[i].b.value + prime_field::mod) % prime_field::mod;

			add_mult_sum[i].b.value = (add_mult_sum[g_zero].a.value * previous_random.value + add_mult_sum[g_zero].b.value) % prime_field::mod;
			add_mult_sum[i].a.value = (add_mult_sum[g_one].a.value * previous_random.value + add_mult_sum[g_one].b.value - add_mult_sum[i].b.value + prime_field::mod) % prime_field::mod;

		}
		ret.a.value = (ret.a.value + add_mult_sum[i].a.value * V_mult_add[i].a.value) % prime_field::mod;
		ret.b.value = (ret.b.value + add_mult_sum[i].a.value * V_mult_add[i].b.value + add_mult_sum[i].b.value * V_mult_add[i].a.value
									+ addV_array[i].a.value) % prime_field::mod;
		ret.c.value = (ret.c.value + add_mult_sum[i].b.value * V_mult_add[i].b.value
									+ addV_array[i].b.value) % prime_field::mod;
	}

	total_uv >>= 1;
	ret.a.value = ret.a.value % prime_field::mod;
	ret.b.value = ret.b.value % prime_field::mod;
	ret.c.value = ret.c.value % prime_field::mod;
	
	ret.a.value = (ret.a.value + prime_field::mod) % prime_field::mod;
	ret.b.value = (ret.b.value + prime_field::mod) % prime_field::mod;
	ret.c.value = (ret.c.value + prime_field::mod) % prime_field::mod;
	
	
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	total_time += time_span.count();
	return ret;
}

void prover::sumcheck_phase2_init(prime_field::field_element previous_random, const prime_field::field_element* r_u, const prime_field::field_element* one_minus_r_u)
{
	std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
	v_u = V_mult_add[0].eval(previous_random);

	int first_half = length_u >> 1, second_half = length_u - first_half;

	beta_u_fhalf[0] = prime_field::field_element(1);
	beta_u_shalf[0] = prime_field::field_element(1);
	for(int i = 0; i < first_half; ++i)
	{
		for(int j = 0; j < (1 << i); ++j)
		{
			beta_u_fhalf[j | (1 << i)].value = beta_u_fhalf[j].value * r_u[i].value % prime_field::mod;
			beta_u_fhalf[j].value = beta_u_fhalf[j].value * one_minus_r_u[i].value % prime_field::mod;
		}
	}

	for(int i = 0; i < second_half; ++i)
	{
		for(int j = 0; j < (1 << i); ++j)
		{
			beta_u_shalf[j | (1 << i)].value = beta_u_shalf[j].value * r_u[i + first_half].value % prime_field::mod;
			beta_u_shalf[j].value = beta_u_shalf[j].value * one_minus_r_u[i + first_half].value % prime_field::mod;
		}
	}

	int mask_fhalf = (1 << first_half) - 1;
	for(int i = 0; i < (1 << length_u); ++i)
	{
		beta_u[i].value = beta_u_fhalf[i & mask_fhalf].value * beta_u_shalf[i >> first_half].value % prime_field::mod;
	}

	
	total_uv = (1 << C.circuit[sumcheck_layer_id - 1].bit_length);
	int total_g = (1 << C.circuit[sumcheck_layer_id].bit_length);
	prime_field::field_element zero = prime_field::field_element(0);
	for(int i = 0; i < total_uv; ++i)
	{
		add_mult_sum[i].a = zero;
		add_mult_sum[i].b = zero;
		addV_array[i].a = zero;
		addV_array[i].b = zero;
		V_mult_add[i] = circuit_value[sumcheck_layer_id - 1][i];
	}
	int first_g_half = (length_g >> 1);
	int mask_g_fhalf = (1 << (length_g >> 1)) - 1;
	for(int i = 0; i < total_g; ++i)
	{
		int ty = C.circuit[sumcheck_layer_id].gates[i].ty;
		int u = C.circuit[sumcheck_layer_id].gates[i].u;
		int v = C.circuit[sumcheck_layer_id].gates[i].v;
		switch(ty)
		{
			case 1: //mult gate
			{
				auto tmp_u = beta_u_fhalf[u & mask_fhalf].value * beta_u_shalf[u >> first_half].value % prime_field::mod;
				auto tmp_g = (beta_g_r0_fhalf[i & mask_g_fhalf].value * beta_g_r0_shalf[i >> first_g_half].value 
								+ beta_g_r1_fhalf[i & mask_g_fhalf].value * beta_g_r1_shalf[i >> first_g_half].value) % prime_field::mod;
				add_mult_sum[v].b.value = add_mult_sum[v].b.value + (tmp_g * tmp_u % prime_field::mod * v_u.value);
				add_mult_sum[v].b.value = add_mult_sum[v].b.value % prime_field::mod;
				break;
			}
			case 0: //add gate
			{
				auto tmp_u = beta_u_fhalf[u & mask_fhalf].value * beta_u_shalf[u >> first_half].value % prime_field::mod;
				auto tmp_g = (beta_g_r0_fhalf[i & mask_g_fhalf].value * beta_g_r0_shalf[i >> first_g_half].value 
								+ beta_g_r1_fhalf[i & mask_g_fhalf].value * beta_g_r1_shalf[i >> first_g_half].value) % prime_field::mod;
				add_mult_sum[v].b.value = (add_mult_sum[v].b.value + tmp_g * tmp_u);
				addV_array[v].b.value = ((tmp_g * tmp_u % prime_field::mod) * v_u.value + addV_array[v].b.value);

				add_mult_sum[v].b.value = add_mult_sum[v].b.value % prime_field::mod;
				addV_array[v].b.value = addV_array[v].b.value % prime_field::mod;
				break;
			}
			case 5: //sum gate
			{
				auto tmp_g = (beta_g_r0_fhalf[i & mask_g_fhalf].value * beta_g_r0_shalf[i >> first_g_half].value 
							+ beta_g_r1_fhalf[i & mask_g_fhalf].value * beta_g_r1_shalf[i >> first_g_half].value) % prime_field::mod;
				auto tmp_g_vu = tmp_g * v_u.value % prime_field::mod;
				for(int j = u; j < v; ++j)
				{
					auto tmp_u = beta_u_fhalf[j & mask_fhalf].value * beta_u_shalf[j >> first_half].value % prime_field::mod;
					addV_array[0].b.value = (addV_array[0].b.value + tmp_g_vu * tmp_u) % prime_field::mod;
				}
				break;
			}
			case 6: //not gate
			{
				auto tmp_u = beta_u_fhalf[u & mask_fhalf].value * beta_u_shalf[u >> first_half].value % prime_field::mod;
				auto tmp_g = (beta_g_r0_fhalf[i & mask_g_fhalf].value * beta_g_r0_shalf[i >> first_g_half].value 
								+ beta_g_r1_fhalf[i & mask_g_fhalf].value * beta_g_r1_shalf[i >> first_g_half].value) % prime_field::mod;
				addV_array[v].b.value = (addV_array[v].b.value + tmp_g * tmp_u % prime_field::mod + prime_field::mod - tmp_g * tmp_u % prime_field::mod * v_u.value % prime_field::mod) % prime_field::mod;
				break;
			}
			case 7: //minus gate
			{
				auto tmp_u = beta_u_fhalf[u & mask_fhalf].value * beta_u_shalf[u >> first_half].value % prime_field::mod;
				auto tmp_g = (beta_g_r0_fhalf[i & mask_g_fhalf].value * beta_g_r0_shalf[i >> first_g_half].value 
								+ beta_g_r1_fhalf[i & mask_g_fhalf].value * beta_g_r1_shalf[i >> first_g_half].value) % prime_field::mod;
				add_mult_sum[v].b.value = (add_mult_sum[v].b.value + prime_field::mod - tmp_g * tmp_u % prime_field::mod) % prime_field::mod;
				addV_array[v].b.value = ((tmp_g * tmp_u % prime_field::mod) * v_u.value + addV_array[v].b.value) % prime_field::mod;
				break;
			}
			case 8: //xor gate
			{
				auto tmp_u = beta_u_fhalf[u & mask_fhalf].value * beta_u_shalf[u >> first_half].value % prime_field::mod;
				auto tmp_g = (beta_g_r0_fhalf[i & mask_g_fhalf].value * beta_g_r0_shalf[i >> first_g_half].value 
								+ beta_g_r1_fhalf[i & mask_g_fhalf].value * beta_g_r1_shalf[i >> first_g_half].value) % prime_field::mod;
				auto tmp = tmp_g * tmp_u % prime_field::mod;
				add_mult_sum[v].b.value = (add_mult_sum[v].b.value + tmp + prime_field::mod - (prime_field::field_element(2).value * v_u.value) % prime_field::mod * tmp % prime_field::mod) % prime_field::mod;
				addV_array[v].b.value = (addV_array[v].b.value + tmp * v_u.value % prime_field::mod) % prime_field::mod;
				break;
			}
			case 9: //NAAB gate
			{
				auto tmp_u = beta_u_fhalf[u & mask_fhalf].value * beta_u_shalf[u >> first_half].value % prime_field::mod;
				auto tmp_g = (beta_g_r0_fhalf[i & mask_g_fhalf].value * beta_g_r0_shalf[i >> first_g_half].value 
								+ beta_g_r1_fhalf[i & mask_g_fhalf].value * beta_g_r1_shalf[i >> first_g_half].value) % prime_field::mod;
				auto tmp = tmp_g * tmp_u % prime_field::mod;
				add_mult_sum[v].b.value = (add_mult_sum[v].b.value + tmp + prime_field::mod - v_u.value * tmp % prime_field::mod) % prime_field::mod;
				break;
			}
			case 10: //relay gate
			{
				auto tmp_u = beta_u_fhalf[u & mask_fhalf].value * beta_u_shalf[u >> first_half].value % prime_field::mod;
				auto tmp_g = (beta_g_r0_fhalf[i & mask_g_fhalf].value * beta_g_r0_shalf[i >> first_g_half].value 
								+ beta_g_r1_fhalf[i & mask_g_fhalf].value * beta_g_r1_shalf[i >> first_g_half].value) % prime_field::mod;
				auto tmp = tmp_g * tmp_u % prime_field::mod;
				addV_array[v].b.value = (addV_array[v].b.value + tmp * v_u.value) % prime_field::mod;
				assert(v == 0);
				break;
			}
			default:
			{
				break;
			}
		}
	}

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	total_time += time_span.count();
}

quadratic_poly prover::sumcheck_phase2_update(prime_field::field_element previous_random, int current_bit)
{
	std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
	quadratic_poly ret = quadratic_poly(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));
	for(int i = 0; i < (total_uv >> 1); ++i)
	{
		int g_zero = i << 1, g_one = i << 1 | 1;
		if(current_bit == 0)
		{
			V_mult_add[i].b = V_mult_add[g_zero].b;
			V_mult_add[i].a = V_mult_add[g_one].b - V_mult_add[i].b;

			addV_array[i].b = addV_array[g_zero].b;
			addV_array[i].a = addV_array[g_one].b - addV_array[i].b;

			add_mult_sum[i].b = add_mult_sum[g_zero].b;
			add_mult_sum[i].a = add_mult_sum[g_one].b - add_mult_sum[i].b;
		}
		else
		{
			
			V_mult_add[i].b.value = (V_mult_add[g_zero].a.value * previous_random.value + V_mult_add[g_zero].b.value) % prime_field::mod;
			V_mult_add[i].a.value = (V_mult_add[g_one].a.value * previous_random.value + V_mult_add[g_one].b.value + prime_field::mod - V_mult_add[i].b.value) % prime_field::mod;

			addV_array[i].b.value = (addV_array[g_zero].a.value * previous_random.value + addV_array[g_zero].b.value) % prime_field::mod;
			addV_array[i].a.value = (addV_array[g_one].a.value * previous_random.value + addV_array[g_one].b.value + prime_field::mod - addV_array[i].b.value) % prime_field::mod;

			add_mult_sum[i].b.value = (add_mult_sum[g_zero].a.value * previous_random.value + add_mult_sum[g_zero].b.value) % prime_field::mod;
			add_mult_sum[i].a.value = (add_mult_sum[g_one].a.value * previous_random.value + add_mult_sum[g_one].b.value + prime_field::mod - add_mult_sum[i].b.value) % prime_field::mod;
		}

		ret.a.value = (ret.a.value + add_mult_sum[i].a.value * V_mult_add[i].a.value) % prime_field::mod;
		ret.b.value = (ret.b.value + add_mult_sum[i].a.value * V_mult_add[i].b.value
									+ add_mult_sum[i].b.value * V_mult_add[i].a.value
									+ addV_array[i].a.value) % prime_field::mod;
		ret.c.value = (ret.c.value + add_mult_sum[i].b.value * V_mult_add[i].b.value
									+ addV_array[i].b.value) % prime_field::mod;
	}

	total_uv >>= 1;

	ret.a.value = ret.a.value % prime_field::mod;
	ret.b.value = ret.b.value % prime_field::mod;
	ret.c.value = ret.c.value % prime_field::mod;
	ret.a.value = (ret.a.value + prime_field::mod) % prime_field::mod;
	ret.b.value = (ret.b.value + prime_field::mod) % prime_field::mod;
	ret.c.value = (ret.c.value + prime_field::mod) % prime_field::mod;
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	total_time += time_span.count();
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

