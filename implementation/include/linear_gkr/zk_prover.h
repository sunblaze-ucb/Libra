#ifndef __zk_prover
#define __zk_prover

#include "linear_gkr/circuit_fast_track.h"
#include "linear_gkr/prime_field.h"
#include "linear_gkr/polynomial.h"
#include <cstring>
#include <utility>
#include <vector>
#include <chrono>
#include "bn.h"
#include "VPD/vpdR.h"
#include "VPD/vpd_test.h"
#include "VPD/input_vpd.h"

class zk_prover
{
public:
	prime_field::field_element v_u, v_v;
	int total_uv;
	layered_circuit C;
	prime_field::field_element* circuit_value[1000000];

	int sumcheck_layer_id, length_g, length_u, length_v;
	prime_field::field_element alpha, beta;
	const prime_field::field_element *r_0, *r_1;
	prime_field::field_element *one_minus_r_0, *one_minus_r_1;

	linear_poly *addV_array;
	linear_poly *V_mult_add;
	//prime_field::field_element *beta_u;
	prime_field::field_element *beta_g_r0_fhalf, *beta_g_r0_shalf, *beta_g_r1_fhalf, *beta_g_r1_shalf, *beta_u_fhalf, *beta_u_shalf;
	prime_field::field_element *beta_u, *beta_v, *beta_g;
	//prime_field::field_element *beta_g_sum;
	linear_poly *add_mult_sum;



	double total_time;
	void init_array(int);
	void get_circuit(const layered_circuit &from_verifier);
	prime_field::field_element* evaluate();
	void proof_init();
	
	std::vector<bn::Ec1> sumcheck_init(int layer_id, int bit_length_g, int bit_length_u, int bit_length_v, const prime_field::field_element &,
		const prime_field::field_element &, const prime_field::field_element*, const prime_field::field_element*, prime_field::field_element*, prime_field::field_element*);
	void sumcheck_phase1_init();
	void sumcheck_phase2_init(prime_field::field_element, const prime_field::field_element*, const prime_field::field_element*);
	quadratic_poly sumcheck_phase1_update(prime_field::field_element, int);
	quintuple_poly sumcheck_phase1_updatelastbit(prime_field::field_element, int);
	quadratic_poly sumcheck_phase2_update(prime_field::field_element, int);
	quintuple_poly sumcheck_phase2_updatelastbit(prime_field::field_element, int);

	quadratic_poly sumcheck_finalround(prime_field::field_element, int, prime_field::field_element);
	prime_field::field_element V_res(const prime_field::field_element*, const prime_field::field_element*, const prime_field::field_element*, int, int);
	std::pair<prime_field::field_element, prime_field::field_element> sumcheck_finalize(prime_field::field_element);
	void delete_self();
	zk_prover()
	{
		memset(circuit_value, 0, sizeof circuit_value);
	}
	~zk_prover();


	//new zk function 


	prime_field::field_element *maskpoly; 
	std::vector<mpz_class> maskpoly_gmp;
	prime_field::field_element maskpoly_sumc;
	prime_field::field_element maskpoly_sumr;
	prime_field::field_element rho;
	std::vector<bn::Ec1> generate_maskpoly_pre_rho(int, int);
	void generate_maskpoly_after_rho(int, int);
	std::pair<std::vector<bn::Ec1>, std::vector<bn::Ec1> > prove_R(std::vector<mpz_class> R, mpz_class &ans);
	std::pair<std::vector<bn::Ec1>, std::vector<bn::Ec1> > prove_mask(std::vector<mpz_class> R, mpz_class &ans);
	std::pair<std::vector<bn::Ec1>, std::vector<bn::Ec1> > prove_input(std::vector<mpz_class> R, mpz_class &ans, mpz_class Z);
	prime_field::field_element query(prime_field::field_element*, prime_field::field_element*, prime_field::field_element);
	prime_field::field_element queryRg1(prime_field::field_element);
	prime_field::field_element queryRg2(prime_field::field_element);

	prime_field::field_element maskR[6], preR[6];
	mpz_class r_f_R, r_f_mask_poly, r_f_input, r_f_input2;
	prime_field::field_element maskR_sumcu, maskR_sumcv, preZu, preZv, Zu, Zv, preu1, prev1, Iuv, prepreu1, preprev1;
	quadratic_poly Rg1, Rg2, sumRc;
	std::vector<bn::Ec1> generate_maskR(int);

	void sumcheck_maskpoly_init();
	std::vector<mpz_class> input_mpz, maskr_mpz;
	std::vector<prime_field::field_element> maskr;
	std::pair<std::vector<bn::Ec1>, std::vector<bn::Ec1> > keygen_and_commit(int input_bit_length, double &key_gen_time);

};

#endif
