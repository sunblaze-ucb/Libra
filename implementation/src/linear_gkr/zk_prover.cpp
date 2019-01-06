#include "linear_gkr/zk_prover.h"

void zk_prover::get_circuit(const layered_circuit &from_verifier)
{
	C = from_verifier;
}

prime_field::field_element zk_prover::V_res(const prime_field::field_element* one_minus_r_0, const prime_field::field_element* r_0, const prime_field::field_element* output_raw, int r_0_size, int output_size)
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

prime_field::field_element* zk_prover::evaluate()
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
		assert(ty == 3);
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
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	std::cerr << "total evaluation time: " << time_span.count() << " seconds." << std::endl;
	return circuit_value[C.total_depth - 1];
}

//new zk function
void zk_prover::generate_maskpoly(int length, int degree){
	if(maskpoly != NULL) 
		delete[] maskpoly;
	maskpoly = new prime_field::field_element[length * degree + 1];
	for(int i = 0; i < length * degree + 1; i++)
		maskpoly[i] = prime_field::random(); 

	maskpoly_sumc = maskpoly[length * degree];
	for(int i = length * degree; i >= 0; i--)
		maskpoly_sumc.value = maskpoly_sumc.value + maskpoly[i].value;
	for(int i = 1; i < length; i++)
		maskpoly_sumc.value = maskpoly_sumc.value + maskpoly_sumc.value;

	//std::cout << "maskpoly_sumc = " << maskpoly_sumc.value << std::endl;
	maskpoly_sumr = prime_field::field_element(0);
}

//new zk function
void zk_prover::sumcheck_init(int layer_id, int bit_length_g, int bit_length_u, int bit_length_v, 
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
	int sumcheck_length = length_u + length_v;
	std::cout << "sumcheck_length = " << sumcheck_length << std::endl;
	generate_maskpoly(sumcheck_length, 2);
}

//new zk function
void zk_prover::generate_maskR(){
	int depth = this->C.total_depth;
	maskR = new prime_field::field_element*[depth];
	for(int i = 0; i < depth; i++){
		maskR[i] = new prime_field::field_element[9];
		for(int j = 0; j < 9; j++)
			maskR[i][j] = prime_field::random();
	}
}


void zk_prover::init_array(int max_bit_length)
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
	addV_array_counter = new int[(1 << max_bit_length)];
	add_mult_sum_counter = new int[(1 << max_bit_length)]; 
}

void zk_prover::delete_self()
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
	delete[] addV_array_counter;
	delete[] add_mult_sum_counter;
	for(int i = 0; i < C.total_depth; ++i)
		delete[] circuit_value[i];
	delete[] maskpoly;
	for(int i = 0; i < C.total_depth; i++)
		delete[] maskR[i];
	delete[] maskR;
}

zk_prover::~zk_prover()
{
}


void zk_prover::sumcheck_phase1_init()
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
		if(C.circuit[sumcheck_layer_id].gates[i].ty == 0) //add gate
		{
			addV_array[u].b.value = (addV_array[u].b.value + circuit_value[sumcheck_layer_id - 1][v].value * beta_g_sum[i].value) % prime_field::mod;
			add_mult_sum[u].b.value = (add_mult_sum[u].b.value + beta_g_sum[i].value) % prime_field::mod;
		}
		if(C.circuit[sumcheck_layer_id].gates[i].ty == 1) //mult gate
		{
			add_mult_sum[u].b.value = (add_mult_sum[u].b.value + circuit_value[sumcheck_layer_id - 1][v].value * beta_g_sum[i].value) % prime_field::mod;
		}
	}
	
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	double tmp = time_span.count();
	total_time += time_span.count();
}


//new zk function
quadratic_poly zk_prover::sumcheck_phase1_update(prime_field::field_element previous_random, int current_bit)
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
		if(i % 8 == 0 || i + 1 == (total_uv >> 1))
		{
			ret.a.value = (ret.a.value + add_mult_sum[i].a.value * V_mult_add[i].a.value) % prime_field::mod;
			ret.b.value = (ret.b.value + add_mult_sum[i].a.value * V_mult_add[i].b.value + add_mult_sum[i].b.value * V_mult_add[i].a.value
									  + addV_array[i].a.value) % prime_field::mod;
			ret.c.value = (ret.c.value + add_mult_sum[i].b.value * V_mult_add[i].b.value
									  + addV_array[i].b.value) % prime_field::mod;
		}
		else
		{
			ret.a.value = (ret.a.value + add_mult_sum[i].a.value * V_mult_add[i].a.value);
			ret.b.value = (ret.b.value + add_mult_sum[i].a.value * V_mult_add[i].b.value + add_mult_sum[i].b.value * V_mult_add[i].a.value
									  + addV_array[i].a.value);
			ret.c.value = (ret.c.value + add_mult_sum[i].b.value * V_mult_add[i].b.value
									  + addV_array[i].b.value);
		}
	}

	total_uv >>= 1;

	//compute with maskpolyR


 	//to do


	//compute with sumcheck maskpoly
	
	prime_field::field_element tmp1, tmp2;
	tmp1.value = maskpoly[current_bit << 1].value;
	tmp2.value = maskpoly[(current_bit << 1) + 1].value;
	tmp1.value = tmp1.value << (length_u + length_v - current_bit - 1);
	tmp2.value = tmp2.value << (length_u + length_v - current_bit - 1);
	maskpoly_sumc.value = (maskpoly_sumc.value - tmp1.value - tmp2.value) / 2;

	prime_field::field_element tmp3;
	if(current_bit > 0){
		maskpoly_sumr = maskpoly_sumr + maskpoly[(current_bit << 1) - 2] * previous_random * previous_random + maskpoly[(current_bit << 1) - 1] * previous_random; 
		tmp3 = maskpoly_sumr;
		for(int i = 0; i < length_u + length_v - current_bit - 1; i++)
			tmp3 = tmp3 + tmp3;
	}

	ret.a.value = (ret.a.value + tmp1.value)% prime_field::mod;
	ret.b.value = (ret.b.value + tmp2.value)% prime_field::mod;
	ret.c.value = (ret.c.value + maskpoly_sumc.value + tmp3.value)% prime_field::mod;

	
	
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	total_time += time_span.count();
	return ret;
}


void zk_prover::sumcheck_phase2_init(prime_field::field_element previous_random, const prime_field::field_element* r_u, const prime_field::field_element* one_minus_r_u)
{
	std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
	v_u = V_mult_add[0].eval(previous_random);
//	fprintf(stderr, "v_u %s\n", v_u.to_string(10).c_str());

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
		addV_array_counter[i] = 0;
		add_mult_sum_counter[i] = 0;
		V_mult_add[i] = circuit_value[sumcheck_layer_id - 1][i];
	}

	for(int i = 0; i < total_g; ++i)
	{
		int ty = C.circuit[sumcheck_layer_id].gates[i].ty;
		int u = C.circuit[sumcheck_layer_id].gates[i].u;
		int v = C.circuit[sumcheck_layer_id].gates[i].v;
		if(ty == 1) //mult gate
		{
			add_mult_sum[v].b.value = add_mult_sum[v].b.value + (beta_g_sum[i].value * beta_u[u].value % prime_field::mod * v_u.value) % prime_field::mod;
			add_mult_sum_counter[v]++;
			if(add_mult_sum_counter[v] > 30)
			{
				add_mult_sum_counter[v] = 0;
				add_mult_sum[v].b.value = add_mult_sum[v].b.value % prime_field::mod;
			}
		}
		if(ty == 0) //add gate
		{
			add_mult_sum[v].b.value = (add_mult_sum[v].b.value + beta_g_sum[i].value * beta_u[u].value);
			addV_array[v].b.value = ((beta_g_sum[i].value * beta_u[u].value % prime_field::mod) * v_u.value + addV_array[v].b.value);

			add_mult_sum_counter[v]++;
			if(add_mult_sum_counter[v] > 30)
			{
				add_mult_sum_counter[v] = 0;
				add_mult_sum[v].b.value = add_mult_sum[v].b.value % prime_field::mod;
			}

			addV_array_counter[v]++;
			if(addV_array_counter[v] > 30)
			{
				addV_array_counter[v] = 0;
				addV_array[v].b.value = addV_array[v].b.value % prime_field::mod;
			}
		}
	}

	for(int i = 0; i < total_uv; ++i)
	{
		if(add_mult_sum_counter[i])
			add_mult_sum[i].b.value = add_mult_sum[i].b.value % prime_field::mod;
		if(addV_array_counter[i])
			addV_array[i].b.value = addV_array[i].b.value % prime_field::mod;
	}

	//update maskpoly
	maskpoly_sumr = maskpoly_sumr + maskpoly[length_u * 2 - 2] * previous_random * previous_random + maskpoly[length_u * 2 - 1] * previous_random;

//	for(int i = 0; i < total_uv; ++i)
//	{
//		fprintf(stderr, "add %s\n", add_array[i].b.to_string(10).c_str());
//		fprintf(stderr, "addV %s\n", add_array[i].b.to_string(10).c_str());
//		fprintf(stderr, "V_mult_add %s\n", add_array[i].b.to_string(10).c_str());
//	}
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	total_time += time_span.count();
}


//new zk function

quadratic_poly zk_prover::sumcheck_phase2_update(prime_field::field_element previous_random, int current_bit)
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
			V_mult_add[i].a.value = (V_mult_add[g_one].a.value * previous_random.value + V_mult_add[g_one].b.value - V_mult_add[i].b.value) % prime_field::mod;

			addV_array[i].b.value = (addV_array[g_zero].a.value * previous_random.value + addV_array[g_zero].b.value) % prime_field::mod;
			addV_array[i].a.value = (addV_array[g_one].a.value * previous_random.value + addV_array[g_one].b.value - addV_array[i].b.value) % prime_field::mod;

			add_mult_sum[i].b.value = (add_mult_sum[g_zero].a.value * previous_random.value + add_mult_sum[g_zero].b.value) % prime_field::mod;
			add_mult_sum[i].a.value = (add_mult_sum[g_one].a.value * previous_random.value + add_mult_sum[g_one].b.value - add_mult_sum[i].b.value) % prime_field::mod;
		}

		if(i % 8 == 0 || i + 1 == (total_uv >> 1))
		{
			ret.a.value = (ret.a.value + add_mult_sum[i].a.value * V_mult_add[i].a.value) % prime_field::mod;
			ret.b.value = (ret.b.value + add_mult_sum[i].a.value * V_mult_add[i].b.value
									  +	add_mult_sum[i].b.value * V_mult_add[i].a.value
									  + addV_array[i].a.value) % prime_field::mod;
			ret.c.value = (ret.c.value + add_mult_sum[i].b.value * V_mult_add[i].b.value
									  + addV_array[i].b.value) % prime_field::mod;
		}
		else
		{
			ret.a.value = (ret.a.value + add_mult_sum[i].a.value * V_mult_add[i].a.value);
			ret.b.value = (ret.b.value + add_mult_sum[i].a.value * V_mult_add[i].b.value
									  +	add_mult_sum[i].b.value * V_mult_add[i].a.value
									  + addV_array[i].a.value);
			ret.c.value = (ret.c.value + add_mult_sum[i].b.value * V_mult_add[i].b.value
									  + addV_array[i].b.value);
		}
	}

	total_uv >>= 1;

	int current = current_bit + length_u;

	prime_field::field_element tmp1, tmp2;
	tmp1.value = maskpoly[current << 1].value;
	tmp2.value = maskpoly[(current << 1) + 1].value;
	tmp1.value = tmp1.value << (length_u + length_v - current - 1);
	tmp2.value = tmp2.value << (length_u + length_v - current - 1);
	maskpoly_sumc.value = (maskpoly_sumc.value - tmp1.value - tmp2.value) / 2;

	prime_field::field_element tmp3;
	maskpoly_sumr = maskpoly_sumr + maskpoly[(current << 1) - 2] * previous_random * previous_random + maskpoly[(current << 1) - 1] * previous_random; 
	tmp3 = maskpoly_sumr;
	for(int i = 0; i < length_u + length_v - current - 1; i++)
		tmp3 = tmp3 + tmp3;
	

	ret.a.value = (ret.a.value + tmp1.value)% prime_field::mod;
	ret.b.value = (ret.b.value + tmp2.value)% prime_field::mod;
	ret.c.value = (ret.c.value + maskpoly_sumc.value + tmp3.value)% prime_field::mod;

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	total_time += time_span.count();
	return ret;
}


std::pair<prime_field::field_element, prime_field::field_element> zk_prover::sumcheck_finalize(prime_field::field_element previous_random)
{
	v_v = V_mult_add[0].eval(previous_random);
	return std::make_pair(v_u, v_v);
}
void zk_prover::proof_init()
{
	//todo
}

