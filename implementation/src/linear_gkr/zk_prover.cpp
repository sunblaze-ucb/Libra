#include "linear_gkr/zk_prover.h"
#include <iostream>
#include <cstring>
prime_field::field_element from_string(const char* str)
{
	prime_field::field_element ret = prime_field::field_element(0);
	int len = strlen(str);
	for(int i = 0; i < len; ++i)
	{
		int digit = str[i] - '0';
		ret = ret * prime_field::field_element(10) + prime_field::field_element(digit);
	}
	return ret;
}
prime_field::field_element inv_2;
void zk_prover::get_circuit(const layered_circuit &from_verifier)
{
	C = from_verifier;
	inv_2 = from_string("8399054365507916142470402071115866954879789801702376374514189432082785107975");
	//std::cout << "test inv_2 = " << inv_2.to_string() << std::endl;

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

	//cheating prover test
	//ret[0].second =  prime_field::field_element(111);

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	std::cerr << "total evaluation time: " << time_span.count() << " seconds." << std::endl;
	return circuit_value[C.total_depth - 1];
}

//a_0x_0 + a_1x_0^2 + a_2x_1 + a_3x_1^2\codts + a_{2n}
void zk_prover::generate_maskpoly(int length, int degree){
	if(maskpoly != NULL) 
		delete[] maskpoly;
	maskpoly = new prime_field::field_element[length * degree + 1];
	for(int i = 0; i < length * degree + 1; i++){
		maskpoly[i] = prime_field::random(); 
		maskpoly[i].value = maskpoly[i].value % prime_field::mod;
	}
	maskpoly_sumc = maskpoly[length * degree];
	//std::cout << "maskpoly[length * degree] = " << maskpoly[length * degree].value << std::endl;
	for(int i = length * degree; i >= 0; i--)
		maskpoly_sumc = maskpoly_sumc + maskpoly[i];
	for(int i = 1; i < length; i++)
		maskpoly_sumc = maskpoly_sumc + maskpoly_sumc;

	//std::cout << "maskpoly_sumc = " << maskpoly_sumc.value << std::endl;
	maskpoly_sumr = prime_field::field_element(0);
}
//new zk function
//a0 + a1g1 + a2g1^2 + a3c + a4c^2 + a5g1c + a6g1^2c + a7g1c^2 + a8g1^2c^2;
void zk_prover::generate_maskR(int layer_id){
	//summaskR = \sum_{c \in {0, 1}}maskR(g1, u);

	summaskR.a = maskR[2] + maskR[2] + maskR[6] + maskR[8];
	summaskR.b = maskR[1] + maskR[1] + maskR[5] + maskR[7];
	summaskR.c = maskR[0] + maskR[0] + maskR[3] + maskR[4];
	//std::cout << "summaskR.a = " << summaskR.a.to_string() << std::endl;
	maskR_sumcu = alpha * Zu * summaskR.eval(preu1);
	//std::cout << "Zu = " << Zu.to_string() << std::endl;
	//std::cout << "preu1 = " << preu1.to_string() << std::endl;
	//maskR_sumcu.value = (maskR_sumcu.value * alpha.value) % prime_field::mod;
	//std::cout << "maskR_sumcu = " << summaskR.eval(preu1).to_string() << std::endl;
	//std::cout << "alpha = " << alpha.to_string() << std::endl;
	//std::cout << "beta = " << beta.to_string() << std::endl;

	maskR_sumcv = beta * Zv * summaskR.eval(prev1); 
	//std::cout << "Zv = " << Zv.to_string() << std::endl;

	//std::cout << "prev1 = " << prev1.to_string() << std::endl;
	//std::cout << "maskR_sumcv = " << summaskR.eval(prev1).to_string() << std::endl;

	//maskR_sumcv.value = (maskR_sumcv.value * beta.value) % prime_field::mod;
	preZu = Zu;
	preZv = Zv;
	Zu = prime_field::field_element(1);
	Zv = prime_field::field_element(1);
	Iuv = prime_field::field_element(1);
	//quadratic function of c that R(z, c) when z = g1;
	Rg1.a = maskR[4] + maskR[7] * preu1 + maskR[8] * preu1 * preu1;
	Rg1.b = maskR[3] + maskR[5] * preu1 + maskR[6] * preu1 * preu1;
	Rg1.c = maskR[0] + maskR[1] * preu1 + maskR[2] * preu1 * preu1;
	//quadratic function of c that R(z, c) when z = g2;
	//std::cout << "maskR_sumcu = " << (Rg1.a + Rg1.b + Rg1.c + Rg1.c).to_string() << std::endl;

	Rg2.a = maskR[4] + maskR[7] * prev1 + maskR[8] * prev1 * prev1;
	Rg2.b = maskR[3] + maskR[5] * prev1 + maskR[6] * prev1 * prev1;
	Rg2.c = maskR[0] + maskR[1] * prev1 + maskR[2] * prev1 * prev1;
	//std::cout << "maskR_sumcv = " << (Rg2.a + Rg2.b + Rg2.c + Rg2.c).to_string() << std::endl;

	for(int i = 0; i < 9; i++)
		premaskR[i] = maskR[i];
	if(layer_id > 1){
		for(int i = 0; i < 9; i++)
			maskR[i] = prime_field::random();
		sumRc.a = maskR[2] + maskR[2] + maskR[6] + maskR[8];
		sumRc.b = maskR[1] + maskR[1] + maskR[5] + maskR[7];
		sumRc.c = maskR[0] + maskR[0] + maskR[3] + maskR[4];
	} 
	if(layer_id == 1){
		maskR[0] = prime_field::random();
		maskR[1] = prime_field::random();
		sumRc.a = prime_field::field_element(0);
		sumRc.b = maskR[1];
		sumRc.c = maskR[0];
	}
}

prime_field::field_element zk_prover::query(prime_field::field_element *u, prime_field::field_element *v, prime_field::field_element r_c){
	//std::cout << "u[1].value = " << u[1].value << std::endl;
	prime_field::field_element result;
	//std::cout << "length_u = " << length_u << "length_v = " << length_v << std::endl; 
	for(int i = 0; i < length_u; i++)
		result = result + maskpoly[2*i] * u[i] * u[i] + maskpoly[2*i + 1] * u[i];
	if(result.value > prime_field::mod) std::cout << "overflow!!!" << std::endl;
	for(int i = 0; i < length_v; i++)
		result = result + maskpoly[2*(i + length_u)] * v[i] * v[i] + maskpoly[2 * (i + length_u) + 1] * v[i];
	//if(result.value > prime_field::mod) std::cout << "overflow!!!" << std::endl;
	//std::cout << "maskpoly[2 * (length_u + length_v)] = " <<  maskpoly[2 * (length_u + length_v)].value << std::endl;
	//if(maskpoly[2 * (length_u + length_v)].value > prime_field::mod)
	//	std::cout << "maskpoly[2 * (length_u + length_v)] = " <<  maskpoly[2 * (length_u + length_v)].value << std::endl;
	result = result + maskpoly[2 * (length_u + length_v)] * r_c * r_c + maskpoly[2 * (length_u + length_v) + 1] * r_c;
	result = result + maskpoly[2 * (length_u + length_v + 1)];
	//if(result.value > prime_field::mod) std::cout << "here overflow!!!" << result.value << std::endl;

	return result;
}	

prime_field::field_element zk_prover::queryRg1(prime_field::field_element r_c){
	return Iuv * preZu * Rg1.eval(r_c) * alpha;
}

prime_field::field_element zk_prover::queryRg2(prime_field::field_element r_c){
	return Iuv * preZv * Rg2.eval(r_c) * beta;	
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
	//std::cout << "sumcheck_length = " << sumcheck_length << std::endl;
	generate_maskpoly(sumcheck_length + 1, 2);
	generate_maskR(layer_id);
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
	//std::cout << "here error" << std::endl;
	//for(int i = 0; i < C.total_depth; i++)
	//	delete[] maskR[i];
	//delete[] maskR;
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
	//quintuple_poly ret5(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));

	//std::cout << "current_bit = " << current_bit << std::endl;
	//std::cout << "total_uv = " << total_uv << std::endl;

	//if(current_bit == 1) preu1 = previous_random;
	//prime_field::field_element tmp5, tmp6;
	Iuv = Iuv * (prime_field::field_element(1) - previous_random);
	//std::cout << "Iuv = " << Iuv.to_string() << std::endl;
	if(current_bit > 0){
		maskR_sumcu = maskR_sumcu * (prime_field::field_element(1) - previous_random);
		maskR_sumcv = maskR_sumcv * (prime_field::field_element(1) - previous_random);

		Zu = Zu * (prime_field::field_element(1) - previous_random) * previous_random;
	}
	//std::cout << "test value" << std::endl;

	//std::cout << "maskR_sumcu = " << maskR_sumcu.to_string() << std::endl;
	//std::cout << "maskR_sumcv = " << maskR_sumcv.value << std::endl;

	ret.b = ret.b - maskR_sumcu - maskR_sumcv;
	ret.c = ret.c + maskR_sumcu + maskR_sumcv;


	
	//compute with sumcheck maskpol
	//quintuple_poly ret5(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));
	

	prime_field::field_element tmp1, tmp2;
	tmp1.value = maskpoly[current_bit << 1].value;
	tmp2.value = maskpoly[(current_bit << 1) + 1].value;

	for(int i = 0; i < length_u + length_v - current_bit; i++){
		tmp1 = tmp1 + tmp1;
		tmp2 = tmp2 + tmp2;
	}

	//tmp1.value = tmp1.value << (length_u + length_v - current_bit - 1);
	//tmp2.value = tmp2.value << (length_u + length_v - current_bit - 1);
	//std::cout << "test value" << std::endl;

	//std::cout << "inv_2 = " << inv_2.to_string() << std::endl;

	maskpoly_sumc = (maskpoly_sumc - tmp1 - tmp2) * inv_2;

	prime_field::field_element tmp3;
	if(current_bit > 0){
		maskpoly_sumr = maskpoly_sumr + maskpoly[(current_bit << 1) - 2] * previous_random * previous_random + maskpoly[(current_bit << 1) - 1] * previous_random; 
		tmp3 = maskpoly_sumr;
		for(int i = 0; i < length_u + length_v - current_bit; i++)
			tmp3 = tmp3 + tmp3;
	}

	//ret.a.value = (ret.a.value + tmp1.value)% prime_field::mod;
	//ret.b.value = (ret.b.value + tmp2.value)% prime_field::mod;
	//ret.c.value = (ret.c.value + maskpoly_sumc.value + tmp3.value)% prime_field::mod;

	ret.a = ret.a + tmp1;
	ret.b = ret.b + tmp2;
	ret.c = ret.c + maskpoly_sumc + tmp3;
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	total_time += time_span.count();
	return ret;
}

quintuple_poly zk_prover::sumcheck_phase1_updatelastbit(prime_field::field_element previous_random, int current_bit)
{	
	std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
	quintuple_poly ret = quintuple_poly(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));
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
	//quintuple_poly ret5(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));

	//std::cout << "current_bit = " << current_bit << std::endl;
	//std::cout << "total_uv = " << total_uv << std::endl;

	//if(current_bit == 1) preu1 = previous_random;
	//prime_field::field_element tmp5, tmp6;
	Iuv = Iuv * (prime_field::field_element(1) - previous_random);
		//std::cout << "Iuv = " << Iuv.to_string() << std::endl;

	if(current_bit > 0){
		maskR_sumcv = maskR_sumcv * (prime_field::field_element(1) - previous_random);
		maskR_sumcu = maskR_sumcu * (prime_field::field_element(1) - previous_random);
		Zu = Zu * (prime_field::field_element(1) - previous_random) * previous_random;
	}
	//std::cout << "maskR_sumcu = " << maskR_sumcu.to_string() << std::endl;
	//std::cout << "maskR_sumcv = " << maskR_sumcv.value << std::endl;

	ret.b = ret.b - maskR_sumcu - maskR_sumcv;
	ret.c = ret.c + maskR_sumcu + maskR_sumcv;
	
	//compute with sumcheck maskpol
	//quintuple_poly ret5(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));
	
	if(current_bit == length_u - 1){
		//std::cout << "come in" << std::endl;
		//(du + e)u(1 - u)(au^2 + bu + c) 
		//(-du^3 + (d - e)u^2 + eu)(au^2 + bu + c)
		//-adu^5 + (a * (d - e) - db)u^4 + (ae + b(d-e) -cd)u^3 + (c(d - e) + be)u^2 + ceu
		prime_field::field_element a = sumRc.a;
		prime_field::field_element b = sumRc.b;
		prime_field::field_element c = sumRc.c;
		prime_field::field_element d = add_mult_sum[0].a;
		prime_field::field_element e = add_mult_sum[0].b;
		//d: coeff of x^5; e coeff of x^4; f coeff of x^3
		ret.d = (prime_field::field_element(0) - a) * d * Zu;
		ret.e = (a * (d - e) - b * d) * Zu;
		ret.f = (a * e + b * (d - e) - c * d) * Zu;
		ret.a = ret.a + (c * (d - e) + b * e) * Zu;  
		ret.b = ret.b + c * e * Zu;
		//prime_field::field_element tmp = ret.d + ret.e + ret.f + (c * (d - e) + b * e) * Zu + c * e * Zu;
		//std::cout << "tmp = " << tmp.to_string() << std::endl;

		//prime_field::field_element tmp1 = (prime_field::field_element(0) - a) * d;
		//if(d prime_field::mod) std::cout << "ddddd" << std::endl;
		//std::cout << "tmp1 = " << tmp1.to_string() << std::endl;
		//prime_field::field_element tmp2 = a * d;
		//std::cout << "tmp2 = " << tmp2.to_string() << std::endl;
		//prime_field::field_element tmp3 = tmp1 + tmp2;
		//std::cout << "tmp3 = " << tmp3.to_string() << std::endl;
		//std::cout << "ret.d= " << ret.d.to_string() << std::endl;
		//std::cout << "ret.e= " << ret.e.to_string() << std::endl;
		//std::cout << "ret.f= " << ret.f.to_string() << std::endl;

	}
	
	prime_field::field_element tmp1, tmp2;
	tmp1.value = maskpoly[current_bit << 1].value;
	tmp2.value = maskpoly[(current_bit << 1) + 1].value;
	for(int i = 0; i < length_u + length_v - current_bit; i++){
		tmp1 = tmp1 + tmp1;
		tmp2 = tmp2 + tmp2;
	}
	//tmp1.value = tmp1.value << (length_u + length_v - current_bit - 1);
	//tmp2.value = tmp2.value << (length_u + length_v - current_bit - 1);
	maskpoly_sumc = (maskpoly_sumc - tmp1 - tmp2) * inv_2;

	prime_field::field_element tmp3;
	if(current_bit > 0){
		maskpoly_sumr = maskpoly_sumr + maskpoly[(current_bit << 1) - 2] * previous_random * previous_random + maskpoly[(current_bit << 1) - 1] * previous_random; 
		tmp3 = maskpoly_sumr;
		for(int i = 0; i < length_u + length_v - current_bit; i++)
			tmp3 = tmp3 + tmp3;
	}

	ret.a.value = (ret.a.value + tmp1.value)% prime_field::mod;
	ret.b.value = (ret.b.value + tmp2.value)% prime_field::mod;
	ret.c.value = (ret.c.value + maskpoly_sumc.value + tmp3.value)% prime_field::mod;

	quintuple_poly ret5(ret.d, ret.e, ret.f, ret.a, ret.b, ret.c);
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	total_time += time_span.count();
	return ret5;
}

void zk_prover::sumcheck_phase2_init(prime_field::field_element previous_random, const prime_field::field_element* r_u, const prime_field::field_element* one_minus_r_u)
{
	std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
	preu1 = previous_random;
	maskR_sumcu = maskR_sumcu * (prime_field::field_element(1) - previous_random);
	maskR_sumcv = maskR_sumcv * (prime_field::field_element(1) - previous_random);

	Iuv = Iuv * (prime_field::field_element(1) - previous_random);
	v_u = V_mult_add[0].eval(previous_random);
	//update v_u to \dot(v_u);
	Zu = Zu * (prime_field::field_element(1) - previous_random) * previous_random; 
	v_u = v_u + Zu * sumRc.eval(previous_random);
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
	//maskR
	//std::cout << "current_bit = " << current_bit << std::endl;
	//if(current_bit == 1) preu1 = previous_random;
	//prime_field::field_element tmp5, tmp6;
	if(current_bit > 0) 
		Iuv = Iuv * (prime_field::field_element(1) - previous_random);
	//std::cout << "Iuv = " << Iuv.to_string() << std::endl;

	if(current_bit > 0){
		//Iuv = Iuv * (prime_field::field_element(1) - previous_random);
		maskR_sumcu = maskR_sumcu * (prime_field::field_element(1) - previous_random);
		maskR_sumcv = maskR_sumcv * (prime_field::field_element(1) - previous_random);
		Zv = Zv * (prime_field::field_element(1) - previous_random) * previous_random;
	}
	//std::cout << "maskR_sumcv = " << maskR_sumcv.to_string() << std::endl;
	ret.b = ret.b - maskR_sumcu - maskR_sumcv;
	ret.c = ret.c + maskR_sumcu + maskR_sumcv;


	//mask sumcheck
	int current = current_bit + length_u;

	prime_field::field_element tmp1, tmp2;
	tmp1.value = maskpoly[current << 1].value;
	tmp2.value = maskpoly[(current << 1) + 1].value;
	for(int i = 0; i < length_u + length_v - current; i++){
		tmp1 = tmp1 + tmp1;
		tmp2 = tmp2 + tmp2;
	}
	//tmp1.value = tmp1.value << (length_u + length_v - current - 1);
	//tmp2.value = tmp2.value << (length_u + length_v - current - 1);
	maskpoly_sumc = (maskpoly_sumc - tmp1 - tmp2) * inv_2;

	prime_field::field_element tmp3;
	maskpoly_sumr = maskpoly_sumr + maskpoly[(current << 1) - 2] * previous_random * previous_random + maskpoly[(current << 1) - 1] * previous_random; 
	tmp3 = maskpoly_sumr;
	for(int i = 0; i < length_u + length_v - current; i++)
		tmp3 = tmp3 + tmp3;
	

	ret.a.value = (ret.a.value + tmp1.value) % prime_field::mod;
	ret.b.value = (ret.b.value + tmp2.value) % prime_field::mod;
	ret.c.value = (ret.c.value + maskpoly_sumc.value + tmp3.value) % prime_field::mod;
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	total_time += time_span.count();
	return ret;
}

quintuple_poly zk_prover::sumcheck_phase2_updatelastbit(prime_field::field_element previous_random, int current_bit)
{
	std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
	quintuple_poly ret = quintuple_poly(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));
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
	//maskR
	//std::cout << "current_bit = " << current_bit << std::endl;
	//if(current_bit == 1) preu1 = previous_random;
	//prime_field::field_element tmp5, tmp6;
	Iuv = Iuv * (prime_field::field_element(1) - previous_random);
	//std::cout << "Iuv = " << Iuv.to_string() << std::endl;

	if(current_bit > 0){
		//Iuv = Iuv * (prime_field::field_element(1) - previous_random);
		maskR_sumcu = maskR_sumcu * (prime_field::field_element(1) - previous_random);

		maskR_sumcv = maskR_sumcv * (prime_field::field_element(1) - previous_random);
		Zv = Zv * (prime_field::field_element(1) - previous_random) * previous_random;
	}
	ret.b = ret.b - maskR_sumcu - maskR_sumcv;
	ret.c = ret.c + maskR_sumcu + maskR_sumcv;
	//std::cout << "maskR_sumcv = " << maskR_sumcv.to_string() << std::endl;
	if(current_bit == length_v - 1){
		//(du + e)u(1 - u)(au^2 + bu + c) 
		//(-du^3 + (d - e)u^2 + eu)(au^2 + bu + c)
		//-adu^5 + (a * (d - e) - db)u^4 + (ae + b(d-e) -cd)u^3 + (c(d - e) + be)u^2 + ceu
		prime_field::field_element a = sumRc.a;
		prime_field::field_element b = sumRc.b;
		prime_field::field_element c = sumRc.c;
		prime_field::field_element d = add_mult_sum[0].a;
		prime_field::field_element e = add_mult_sum[0].b;
		//d: coeff of x^5; e coeff of x^4; f coeff of x^3
		ret.d = (prime_field::field_element(0) - a) * d * Zv;
		ret.e = (a * (d - e) - b * d) * Zv;
		ret.f = (a * e + b * (d - e) - c * d) * Zv;
		ret.a = ret.a + (c * (d - e) + b * e) * Zv;  
		ret.b = ret.b + c * e * Zv;
	}
	//mask sumcheck
	int current = current_bit + length_u;

	prime_field::field_element tmp1, tmp2;
	tmp1.value = maskpoly[current << 1].value;
	tmp2.value = maskpoly[(current << 1) + 1].value;
	for(int i = 0; i < length_u + length_v - current; i++){
		tmp1 = tmp1 + tmp1;
		tmp2 = tmp2 + tmp2;
	}
	//tmp1.value = tmp1.value << (length_u + length_v - current - 1);
	//tmp2.value = tmp2.value << (length_u + length_v - current - 1);
	maskpoly_sumc = (maskpoly_sumc - tmp1 - tmp2) * inv_2;

	prime_field::field_element tmp3;
	maskpoly_sumr = maskpoly_sumr + maskpoly[(current << 1) - 2] * previous_random * previous_random + maskpoly[(current << 1) - 1] * previous_random; 
	tmp3 = maskpoly_sumr;
	for(int i = 0; i < length_u + length_v - current; i++)
		tmp3 = tmp3 + tmp3;
	

	ret.a.value = (ret.a.value + tmp1.value) % prime_field::mod;
	ret.b.value = (ret.b.value + tmp2.value) % prime_field::mod;
	ret.c.value = (ret.c.value + maskpoly_sumc.value + tmp3.value) % prime_field::mod;
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	quintuple_poly ret5(ret.d, ret.e, ret.f, ret.a, ret.b, ret.c);
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	total_time += time_span.count();
	return ret5;
}

quadratic_poly zk_prover::sumcheck_finalround(prime_field::field_element previous_random, int current, prime_field::field_element general_value){
	//prev1 = previous_random;
	quadratic_poly ret = quadratic_poly(prime_field::field_element(0), prime_field::field_element(0), prime_field::field_element(0));
	//to do
	ret.a = Iuv * preZu * Rg1.a * alpha + Iuv * preZv * Rg2.a * beta;
	//std::cout << "ret.a = " << ret.a.to_string() << std::endl;
	ret.b = Iuv * preZu * Rg1.b * alpha + Iuv * preZv * Rg2.b * beta + general_value;
	//std::cout << "alpha = " << alpha.to_string() << std::endl;
	//std::cout << "beta = " << beta.to_string() << std::endl;

	//std::cout << "ret.b = " << ret.b.to_string() << std::endl;
	//std::cout << "general_value = " << general_value.to_string() << std::endl;
	//std::cout << "Iuv = " << Iuv.to_string() << std::endl;
	ret.c = Iuv * preZu * Rg1.c * alpha + Iuv * preZv * Rg2.c * beta;
	//std::cout << "ret.c = " << ret.c.to_string() << std::endl;

	//std::cout << "current = " << current << std::endl;
	prime_field::field_element tmp1, tmp2;
	tmp1.value = maskpoly[current << 1].value;
	tmp2.value = maskpoly[(current << 1) + 1].value;
	for(int i = 0; i < length_u + length_v - current; i++){
		tmp1 = tmp1 + tmp1;
		tmp2 = tmp2 + tmp2;
	}
	//tmp1.value = tmp1.value << (length_u + length_v - current - 1);
	//tmp2.value = tmp2.value << (length_u + length_v - current - 1);
	maskpoly_sumc = (maskpoly_sumc - tmp1 - tmp2) * inv_2;

	prime_field::field_element tmp3;
	maskpoly_sumr = maskpoly_sumr + maskpoly[(current << 1) - 2] * previous_random * previous_random + maskpoly[(current << 1) - 1] * previous_random; 
	tmp3 = maskpoly_sumr;
	for(int i = 0; i < length_u + length_v - current; i++)
		tmp3 = tmp3 + tmp3;
	//std::cout << "sumc = " << maskpoly_sumc.to_string() << std::endl;
	//std::cout << "sumc = " << maskpoly[2 * current + 2].to_string() << std::endl;


	ret.a.value = (ret.a.value + tmp1.value) % prime_field::mod;
	//std::cout << "a = " << ret.a.to_string() << std::endl;
	ret.b.value = (ret.b.value + tmp2.value) % prime_field::mod;
	ret.c.value = (ret.c.value + maskpoly_sumc.value + tmp3.value) % prime_field::mod;
	//std::cout << "previous_random = " << previous_random.to_string() << std::endl;
	return ret;
};


std::pair<prime_field::field_element, prime_field::field_element> zk_prover::sumcheck_finalize(prime_field::field_element previous_random)
{
	prev1 = previous_random;
	//std::cout << "previous_random = " << previous_random.to_string() << std::endl;
	Iuv = Iuv * (prime_field::field_element(1) - previous_random);
	v_v = V_mult_add[0].eval(previous_random);
	Zv = Zv * (prime_field::field_element(1) - previous_random) * previous_random;
	//std::cout << "Zv = " << Zv.to_string() << std::endl;

	v_v = v_v + Zv * sumRc.eval(previous_random);
	//std::cout << "compare = " << Zv * sum
	return std::make_pair(v_u, v_v);
}
void zk_prover::proof_init()
{
	//todo
}
