#include "linear_gkr/prover_fast_para_track.h"
#include <iostream>
#define NUM_THREADS 4

void prover::get_circuit(const layered_circuit &from_verifier)
{
	C = from_verifier;
}

void V_res_func(int i, int id, prime_field::field_element *output, prime_field::field_element *noutput, int upper_bound, const prime_field::field_element* one_minus_r_0, const prime_field::field_element *r_0)
{
	int bs = upper_bound / NUM_THREADS;
	int st = id * bs, ed = (id + 1) * bs;
	for(int j = st; j < ed; ++j)
	{
		noutput[j].value = (output[j << 1].value * one_minus_r_0[i].value + output[j << 1 | 1].value * r_0[i].value) % prime_field::mod;
	}
}

prime_field::field_element prover::V_res(const prime_field::field_element* one_minus_r_0, const prime_field::field_element* r_0, const prime_field::field_element* output_raw, int r_0_size, int output_size)
{
	std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
	prime_field::field_element *output, *noutput;
	output = new prime_field::field_element[output_size];
	noutput = new prime_field::field_element[output_size];
	for(int i = 0; i < output_size; ++i)
		output[i] = output_raw[i];
	for(int i = 0; i < r_0_size; ++i)
	{
		if((output_size >> 1) >= NUM_THREADS)
		{
			std::thread threads[NUM_THREADS];
			for(int j = 0; j < NUM_THREADS; ++j)
				threads[j] = std::thread(V_res_func, i, j, output, noutput, output_size >> 1, one_minus_r_0, r_0);
			for(int j = 0; j < NUM_THREADS; ++j)
				threads[j].join();
		}
		else
		{
			for(int j = 0; j < (output_size >> 1); ++j)
			{
				noutput[j].value = (output[j << 1].value * one_minus_r_0[i].value + output[j << 1 | 1].value * r_0[i].value) % prime_field::mod;
			}
		}
		std::swap(output, noutput);
		output_size >>= 1;
	}
	
	prime_field::field_element res = output[0];
	delete[] output;
	delete[] noutput;
	if(res.value < 0)
		res.value = res.value + prime_field::mod;
		
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	total_time += time_span.count();
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
	nadd_mult_sum = new linear_poly[(1 << (max_bit_length))];
	nV_mult_add = new linear_poly[(1 << (max_bit_length))];
	naddV_array = new linear_poly[(1 << (max_bit_length))];
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

void prover::delete_self()
{
	delete[] add_mult_sum;
	delete[] V_mult_add;
	delete[] addV_array;
	delete[] nadd_mult_sum;
	delete[] nV_mult_add;
	delete[] naddV_array;
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
}

prover::~prover()
{
}

void phase1_init_func(prover* p, int id)
{
	int total_uv = (1 << (p -> C.circuit[p -> sumcheck_layer_id - 1].bit_length));
	int first_half = p -> length_g >> 1, second_half = p -> length_g - first_half;
	prime_field::field_element zero = prime_field::field_element(0);
	int layer_id = p -> sumcheck_layer_id;
	
	int bs = (total_uv) / NUM_THREADS;
	if(bs == 0)
		bs = 1;
	int st = id * bs, ed = (id + 1) * bs;
	if(ed > total_uv)
		ed = total_uv;
	
	for(int i = st; i < ed; ++i)
	{
		p -> V_mult_add[i] = p -> circuit_value[layer_id - 1][i];
		p -> addV_array[i].a = zero;
		p -> addV_array[i].b = zero;
		p -> add_mult_sum[i].a = zero;
		p -> add_mult_sum[i].b = zero;
	}
	

	int mask_fhalf = (1 << first_half) - 1;
	bs = (1 << p -> length_g) / NUM_THREADS;
	if(bs == 0)
		bs = 1;
	st = bs * id, ed = bs * (id + 1);
	if(ed > (1 << p -> length_g))
		ed = (1 << p -> length_g);
	for(int i = st; i < ed; ++i)
	{
		p -> beta_g_sum[i].value = (p -> beta_g_r0_fhalf[i & mask_fhalf].value * p -> beta_g_r0_shalf[i >> first_half].value 
							 + p -> beta_g_r1_fhalf[i & mask_fhalf].value * p -> beta_g_r1_shalf[i >> first_half].value) % prime_field::mod;
	}
	
	/*

	for(int i = st; i < ed; ++i)
	{
		int u, v;
		u = p -> C.circuit[p -> sumcheck_layer_id].gates[i].u;
		v = p -> C.circuit[p -> sumcheck_layer_id].gates[i].v;
		//race condition bug here
		//plan a change next day
		if(p -> C.circuit[p -> sumcheck_layer_id].gates[i].ty == 0) //add gate
		{
			p -> addV_array[u].b.value = (p -> addV_array[u].b.value + p -> circuit_value[p -> sumcheck_layer_id - 1][v].value * p -> beta_g_sum[i].value) % prime_field::mod;
			p -> add_mult_sum[u].b.value = (p -> add_mult_sum[u].b.value + p -> beta_g_sum[i].value) % prime_field::mod;
		}
		if(p -> C.circuit[p -> sumcheck_layer_id].gates[i].ty == 1) //mult gate
		{
			p -> add_mult_sum[u].b.value = (p -> add_mult_sum[u].b.value + p -> circuit_value[p -> sumcheck_layer_id - 1][v].value * p -> beta_g_sum[i].value) % prime_field::mod;
		}
	}
	*/
	
	bs = (total_uv) / NUM_THREADS;
	if(bs == 0)
		bs = 1;
	st = id * bs, ed = (id + 1) * bs;
	for(int i = st; i < ed; ++i)
	{
		if(p -> C.circuit[p -> sumcheck_layer_id].u_gates.find(i) != p -> C.circuit[p -> sumcheck_layer_id].u_gates.end())
		{
			int sz = p -> C.circuit[p -> sumcheck_layer_id].u_gates[i].size();
			for(int j = 0; j < sz; ++j)
			{
				int u = i;
				int g = p -> C.circuit[p -> sumcheck_layer_id].u_gates[i][j].second.first;
				int v = p -> C.circuit[p -> sumcheck_layer_id].u_gates[i][j].second.second;
				int ty = p -> C.circuit[p -> sumcheck_layer_id].u_gates[i][j].first;
				if(ty == 0)
				{
					p -> addV_array[u].b.value = (p -> addV_array[u].b.value + p -> circuit_value[p -> sumcheck_layer_id - 1][v].value * p -> beta_g_sum[g].value) % prime_field::mod;
					p -> add_mult_sum[u].b.value = (p -> add_mult_sum[u].b.value + p -> beta_g_sum[g].value) % prime_field::mod;
				}
				if(ty == 1)
				{
					p -> add_mult_sum[u].b.value = (p -> add_mult_sum[u].b.value + p -> circuit_value[p -> sumcheck_layer_id - 1][v].value * p -> beta_g_sum[g].value) % prime_field::mod;
				}
			}
		}
	}
}

void prover::sumcheck_phase1_init()
{
	std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();

	beta_g_r0_fhalf[0] = alpha;
	beta_g_r1_fhalf[0] = beta;
	beta_g_r0_shalf[0] = prime_field::field_element(1);
	beta_g_r1_shalf[0] = prime_field::field_element(1);

	total_uv = (1 << C.circuit[sumcheck_layer_id - 1].bit_length);

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

	std::thread threads[NUM_THREADS];

	for(int i = 0; i < NUM_THREADS; ++i)
		threads[i] = std::thread(phase1_init_func, this, i);
	for(int i = 0; i < NUM_THREADS; ++i)
		threads[i].join();

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	double tmp = time_span.count();
	total_time += time_span.count();
}

void phase1_update_func(prover *p, int id, int total_uv, int current_bit, prime_field::field_element previous_random)
{
	int bs = (total_uv >> 1) / NUM_THREADS;
	if(bs == 0)
		bs = 1;
	int st = id * bs;
	int ed = (id + 1) * bs;
	if(ed > (total_uv >> 1))
	{
		ed = (total_uv >> 1);
	}
	p -> rets[id].a = prime_field::field_element(0);
	p -> rets[id].b = prime_field::field_element(0);
	p -> rets[id].c = prime_field::field_element(0);
	for(int i = st; i < ed; ++i)
	{
		prime_field::field_element zero_value, one_value;
		int g_zero = i << 1, g_one = i << 1 | 1;
		if(current_bit == 0)
		{
			p -> nV_mult_add[i].b = p -> V_mult_add[g_zero].b;
			p -> nV_mult_add[i].a = p -> V_mult_add[g_one].b - p -> nV_mult_add[i].b;

			p -> naddV_array[i].b = p -> addV_array[g_zero].b;
			p -> naddV_array[i].a = p -> addV_array[g_one].b - p -> naddV_array[i].b;

			p -> nadd_mult_sum[i].b = p -> add_mult_sum[g_zero].b;
			p -> nadd_mult_sum[i].a = p -> add_mult_sum[g_one].b - p -> nadd_mult_sum[i].b;

		}
		else
		{
			p -> nV_mult_add[i].b.value = (p -> V_mult_add[g_zero].a.value * previous_random.value + p -> V_mult_add[g_zero].b.value) % prime_field::mod;
			p -> nV_mult_add[i].a.value = (p -> V_mult_add[g_one].a.value * previous_random.value + p -> V_mult_add[g_one].b.value - p -> nV_mult_add[i].b.value + prime_field::mod) % prime_field::mod;

			p -> naddV_array[i].b.value = (p -> addV_array[g_zero].a.value * previous_random.value + p -> addV_array[g_zero].b.value) % prime_field::mod;
			p -> naddV_array[i].a.value = (p -> addV_array[g_one].a.value * previous_random.value + p -> addV_array[g_one].b.value - p -> naddV_array[i].b.value + prime_field::mod) % prime_field::mod;

			p -> nadd_mult_sum[i].b.value = (p -> add_mult_sum[g_zero].a.value * previous_random.value + p -> add_mult_sum[g_zero].b.value) % prime_field::mod;
			p -> nadd_mult_sum[i].a.value = (p -> add_mult_sum[g_one].a.value * previous_random.value + p -> add_mult_sum[g_one].b.value - p -> nadd_mult_sum[i].b.value + prime_field::mod) % prime_field::mod;

		}
		if(i % 4 == 0 || i + 1 == (total_uv >> 1))
		{
			p -> rets[id].a.value = (p -> rets[id].a.value + p -> nadd_mult_sum[i].a.value * p -> nV_mult_add[i].a.value) % prime_field::mod;
			p -> rets[id].b.value = (p -> rets[id].b.value + p -> nadd_mult_sum[i].a.value * p -> nV_mult_add[i].b.value + p -> nadd_mult_sum[i].b.value * p -> nV_mult_add[i].a.value
									  + p -> naddV_array[i].a.value) % prime_field::mod;
			p -> rets[id].c.value = (p -> rets[id].c.value + p -> nadd_mult_sum[i].b.value * p -> nV_mult_add[i].b.value
									  + p -> naddV_array[i].b.value) % prime_field::mod;
		}
		else
		{
			p -> rets[id].a.value = (p -> rets[id].a.value + p -> nadd_mult_sum[i].a.value * p -> nV_mult_add[i].a.value);
			p -> rets[id].b.value = (p -> rets[id].b.value + p -> nadd_mult_sum[i].a.value * p -> nV_mult_add[i].b.value + p -> nadd_mult_sum[i].b.value * p -> nV_mult_add[i].a.value
									  + p -> naddV_array[i].a.value);
			p -> rets[id].c.value = (p -> rets[id].c.value + p -> nadd_mult_sum[i].b.value * p -> nV_mult_add[i].b.value
									  + p -> naddV_array[i].b.value);
		}
	}

	p -> rets[id].a.value = p -> rets[id].a.value % prime_field::mod;
	p -> rets[id].b.value = p -> rets[id].b.value % prime_field::mod;
	p -> rets[id].c.value = p -> rets[id].c.value % prime_field::mod;
}

quadratic_poly prover::sumcheck_phase1_update(prime_field::field_element previous_random, int current_bit)
{
//	fprintf(stderr, "sumcheck level %d, phase1 update start\n", sumcheck_layer_id);
	std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
	quadratic_poly ret = quadratic_poly(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));
	std::thread threads[NUM_THREADS];
	for(int i = 0; i < NUM_THREADS; ++i)
		threads[i] = std::thread(phase1_update_func, this, i, total_uv, current_bit, previous_random);
	for(int i = 0; i < NUM_THREADS; ++i)
		threads[i].join();
	for(int i = 0; i < NUM_THREADS; ++i)
		ret = ret + rets[i];
	
	std::swap(naddV_array, addV_array);
	std::swap(nadd_mult_sum, add_mult_sum);
	std::swap(nV_mult_add, V_mult_add);
	
	total_uv >>= 1;
//	fprintf(stderr, "sumcheck level %d, phase1 update finished\n", sumcheck_layer_id);

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	total_time += time_span.count();
//	std::cerr << "total phase1 update " << current_bit << " time: " << time_span.count() << " seconds" << std::endl;
	return ret;
}

void phase2_init_func(prover* p, int id)
{
	int first_half = p -> length_u >> 1, second_half = p -> length_u - first_half;
	int mask_fhalf = (1 << first_half) - 1;
	
	int bs, st, ed;
	
	bs = (1 << p -> length_u) / NUM_THREADS;
	if(bs == 0)
		bs = 1;
	st = id * bs;
	ed = (id + 1) * bs;
	if(ed > (1 << p -> length_u))
		ed = (1 << p -> length_u);
	
	for(int i = st; i < ed; ++i)
	{
		p -> beta_u[i].value = p -> beta_u_fhalf[i & mask_fhalf].value * p -> beta_u_shalf[i >> first_half].value % prime_field::mod;
	}

	
	int total_uv = (1 << p -> C.circuit[p -> sumcheck_layer_id - 1].bit_length);
	int total_g = (1 << p -> C.circuit[p -> sumcheck_layer_id].bit_length);
	prime_field::field_element zero = prime_field::field_element(0);
	
	bs = total_uv / NUM_THREADS;
	if(bs == 0)
		bs = 1;
	st = id * bs, ed = (id + 1) * bs;
	
	if(ed > total_uv)
		ed = total_uv;
	
	for(int i = st; i < ed; ++i)
	{
		p -> add_mult_sum[i].a = zero;
		p -> add_mult_sum[i].b = zero;
		p -> addV_array[i].a = zero;
		p -> addV_array[i].b = zero;
		p -> addV_array_counter[i] = 0;
		p -> add_mult_sum_counter[i] = 0;
		p -> V_mult_add[i] = p -> circuit_value[p -> sumcheck_layer_id - 1][i];
	}
	
	/*
	bs = total_g / 16;
	if(bs == 0)
		bs = 1;
	st = id * bs, ed = (id + 1) * bs;
	if(ed > total_g)
		ed = total_g;
	for(int i = st; i < ed; ++i)
	{
		int ty = p -> C.circuit[p -> sumcheck_layer_id].gates[i].ty;
		int u = p -> C.circuit[p -> sumcheck_layer_id].gates[i].u;
		int v = p -> C.circuit[p -> sumcheck_layer_id].gates[i].v;
		if(ty == 1) //mult gate
		{
			p -> add_mult_sum[v].b.value = p -> add_mult_sum[v].b.value + (p -> beta_g_sum[i].value * p -> beta_u[u].value % prime_field::mod * p -> v_u.value) % prime_field::mod;
			p -> add_mult_sum_counter[v]++;
			if(p -> add_mult_sum_counter[v] > 30)
			{
				p -> add_mult_sum_counter[v] = 0;
				p -> add_mult_sum[v].b.value = p -> add_mult_sum[v].b.value % prime_field::mod;
			}
		}
		if(ty == 0) //add gate
		{
			p -> add_mult_sum[v].b.value = (p -> add_mult_sum[v].b.value + p -> beta_g_sum[i].value * p -> beta_u[u].value);
			p -> addV_array[v].b.value = ((p -> beta_g_sum[i].value * p -> beta_u[u].value % prime_field::mod) * p -> v_u.value + p -> addV_array[v].b.value);

			p -> add_mult_sum_counter[v]++;
			if(p -> add_mult_sum_counter[v] > 30)
			{
				p -> add_mult_sum_counter[v] = 0;
				p -> add_mult_sum[v].b.value = p -> add_mult_sum[v].b.value % prime_field::mod;
			}

			p -> addV_array_counter[v]++;
			if(p -> addV_array_counter[v] > 30)
			{
				p -> addV_array_counter[v] = 0;
				p -> addV_array[v].b.value = p -> addV_array[v].b.value % prime_field::mod;
			}
		}
	}
	*/
	bs = total_uv / NUM_THREADS;
	if(bs == 0)
		bs = 1;
	st = id * bs;
	ed = (id + 1) * bs;
	if(ed > total_uv)
		ed = total_uv;
	
	for(int i = st; i < ed; ++i)
	{
		if(p -> C.circuit[p -> sumcheck_layer_id].v_gates.find(i) != p -> C.circuit[p -> sumcheck_layer_id].v_gates.end())
		{
			int sz = p -> C.circuit[p -> sumcheck_layer_id].v_gates[i].size();
			for(int j = 0; j < sz; ++j)
			{
				int ty = p -> C.circuit[p -> sumcheck_layer_id].v_gates[i][j].first;
				int g = p -> C.circuit[p -> sumcheck_layer_id].v_gates[i][j].second.first;
				int u = p -> C.circuit[p -> sumcheck_layer_id].v_gates[i][j].second.second;
				int v = i;
				if(ty == 1) //mult gate
				{
					p -> add_mult_sum[v].b.value = p -> add_mult_sum[v].b.value + (p -> beta_g_sum[g].value * p -> beta_u[u].value % prime_field::mod * p -> v_u.value) % prime_field::mod;
					p -> add_mult_sum_counter[v]++;
					if(p -> add_mult_sum_counter[v] > 30)
					{
						p -> add_mult_sum_counter[v] = 0;
						p -> add_mult_sum[v].b.value = p -> add_mult_sum[v].b.value % prime_field::mod;
					}
				}
				if(ty == 0) //add gate
				{
					p -> add_mult_sum[v].b.value = (p -> add_mult_sum[v].b.value + p -> beta_g_sum[g].value * p -> beta_u[u].value);
					p -> addV_array[v].b.value = ((p -> beta_g_sum[g].value * p -> beta_u[u].value % prime_field::mod) * p -> v_u.value + p -> addV_array[v].b.value);

					p -> add_mult_sum_counter[v]++;
					if(p -> add_mult_sum_counter[v] > 30)
					{
						p -> add_mult_sum_counter[v] = 0;
						p -> add_mult_sum[v].b.value = p -> add_mult_sum[v].b.value % prime_field::mod;
					}

					p -> addV_array_counter[v]++;
					if(p -> addV_array_counter[v] > 30)
					{
						p -> addV_array_counter[v] = 0;
						p -> addV_array[v].b.value = p -> addV_array[v].b.value % prime_field::mod;
					}
				}
			}
		}
	}
	
	for(int i = st; i < ed; ++i)
	{
		if(p -> add_mult_sum_counter[i])
			p -> add_mult_sum[i].b.value = p -> add_mult_sum[i].b.value % prime_field::mod;
		if(p -> addV_array_counter[i])
			p -> addV_array[i].b.value = p -> addV_array[i].b.value % prime_field::mod;
	}
}

void prover::sumcheck_phase2_init(prime_field::field_element previous_random, const prime_field::field_element* r_u, const prime_field::field_element* one_minus_r_u)
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
	
	total_uv = (1 << C.circuit[sumcheck_layer_id - 1].bit_length);
	std::thread threads[NUM_THREADS];
	for(int i = 0; i < NUM_THREADS; ++i)
		threads[i] = std::thread(phase2_init_func, this, i);
	
	for(int i = 0; i < NUM_THREADS; ++i)
		threads[i].join();
	
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	total_time += time_span.count();
}

void phase2_update_func(prover *p, int id, int total_uv, int current_bit, prime_field::field_element previous_random)
{
	int bs = (total_uv >> 1) / NUM_THREADS;
	if(bs == 0)
		bs = 1;
	int st = bs * id;
	int ed = bs * (id + 1);
	if(ed > (total_uv >> 1))
		ed = (total_uv >> 1);
	
	for(int i = st; i < ed; ++i)
	{
		int g_zero = i << 1, g_one = i << 1 | 1;
		if(current_bit == 0)
		{
			p -> nV_mult_add[i].b = p -> V_mult_add[g_zero].b;
			p -> nV_mult_add[i].a = p -> V_mult_add[g_one].b - p -> nV_mult_add[i].b;

			p -> naddV_array[i].b = p -> addV_array[g_zero].b;
			p -> naddV_array[i].a = p -> addV_array[g_one].b - p -> naddV_array[i].b;

			p -> nadd_mult_sum[i].b = p -> add_mult_sum[g_zero].b;
			p -> nadd_mult_sum[i].a = p -> add_mult_sum[g_one].b - p -> nadd_mult_sum[i].b;
		}
		else
		{
			
			p -> nV_mult_add[i].b.value = (p -> V_mult_add[g_zero].a.value * previous_random.value + p -> V_mult_add[g_zero].b.value) % prime_field::mod;
			p -> nV_mult_add[i].a.value = (p -> V_mult_add[g_one].a.value * previous_random.value + p -> V_mult_add[g_one].b.value - p -> nV_mult_add[i].b.value) % prime_field::mod;

			p -> naddV_array[i].b.value = (p -> addV_array[g_zero].a.value * previous_random.value + p -> addV_array[g_zero].b.value) % prime_field::mod;
			p -> naddV_array[i].a.value = (p -> addV_array[g_one].a.value * previous_random.value + p -> addV_array[g_one].b.value - p -> naddV_array[i].b.value) % prime_field::mod;

			p -> nadd_mult_sum[i].b.value = (p -> add_mult_sum[g_zero].a.value * previous_random.value + p -> add_mult_sum[g_zero].b.value) % prime_field::mod;
			p -> nadd_mult_sum[i].a.value = (p -> add_mult_sum[g_one].a.value * previous_random.value + p -> add_mult_sum[g_one].b.value - p -> nadd_mult_sum[i].b.value) % prime_field::mod;
		}

		if(i % 8 == 0 || i + 1 == (total_uv >> 1))
		{
			p -> rets[id].a.value = (p -> rets[id].a.value + p -> nadd_mult_sum[i].a.value * p -> nV_mult_add[i].a.value) % prime_field::mod;
			p -> rets[id].b.value = (p -> rets[id].b.value + p -> nadd_mult_sum[i].a.value * p -> nV_mult_add[i].b.value
									  +	p -> nadd_mult_sum[i].b.value * p -> nV_mult_add[i].a.value
									  + p -> naddV_array[i].a.value) % prime_field::mod;
			p -> rets[id].c.value = (p -> rets[id].c.value + p -> nadd_mult_sum[i].b.value * p -> nV_mult_add[i].b.value
									  + p -> naddV_array[i].b.value) % prime_field::mod;
		}
		else
		{
			p -> rets[id].a.value = (p -> rets[id].a.value + p -> nadd_mult_sum[i].a.value * p -> nV_mult_add[i].a.value);
			p -> rets[id].b.value = (p -> rets[id].b.value + p -> nadd_mult_sum[i].a.value * p -> nV_mult_add[i].b.value
									  +	p -> nadd_mult_sum[i].b.value * p -> nV_mult_add[i].a.value
									  + p -> naddV_array[i].a.value);
			p -> rets[id].c.value = (p -> rets[id].c.value + p -> nadd_mult_sum[i].b.value * p -> nV_mult_add[i].b.value
									  + p -> naddV_array[i].b.value);
		}
	}


	p -> rets[id].a.value = p -> rets[id].a.value % prime_field::mod;
	p -> rets[id].b.value = p -> rets[id].b.value % prime_field::mod;
	p -> rets[id].c.value = p -> rets[id].c.value % prime_field::mod;
}

quadratic_poly prover::sumcheck_phase2_update(prime_field::field_element previous_random, int current_bit)
{
//	fprintf(stderr, "sumcheck level %d, phase2 update start\n", sumcheck_layer_id);
	std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
	quadratic_poly ret = quadratic_poly(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));
	
	std::thread threads[NUM_THREADS];
	for(int i = 0; i < NUM_THREADS; ++i)
		rets[i] = quadratic_poly(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));
	for(int i = 0; i < NUM_THREADS; ++i)
		threads[i] = std::thread(phase2_update_func, this, i, total_uv, current_bit, previous_random);
	for(int i = 0; i < NUM_THREADS; ++i)
		threads[i].join();
	
	for(int i = 0; i < NUM_THREADS; ++i)
		ret = ret + rets[i];
		
	std::swap(naddV_array, addV_array);
	std::swap(nadd_mult_sum, add_mult_sum);
	std::swap(nV_mult_add, V_mult_add);
	
	total_uv >>= 1;
//	fprintf(stderr, "sumcheck level %d, phase2 update finished\n", sumcheck_layer_id);
	
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	total_time += time_span.count();
//	std::cerr << "total phase2 update " << current_bit << " time: " << time_span.count() << " seconds" << std::endl;
	
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
