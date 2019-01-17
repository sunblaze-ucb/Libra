#include <cstdlib>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <gmp.h>
#include <gmpxx.h>
#include <math.h>
#include <string>

#include "test_point.hpp"
#include "bn.h"

using namespace std;
using namespace bn;

#define P 512
const int multi_scalar_w = 4;

unsigned long int seed;
int NumOfVar;
gmp_randstate_t r_state;  
mpz_class p, a;
Ec1 g1, g1a, result;
Ec2 g2, g2a;

vector<Ec1> pub_g1, g1_pre;
vector<Ec2> pub_g2, g2_pre;

class multi_scalar_state
{
public:
	Ec1 value[1 << multi_scalar_w];
};

vector<multi_scalar_state> multi_scalar_g1;

std::vector<mpz_class> pub_g1_exp;

//mpz_class globalright, globalleft;

vector<mpz_class> s;

std::vector<mpz_class> pre_input(std::vector<mpz_class>& input){
	cout << "inputsize = " << input.size() << endl; 
	int total = 1 << NumOfVar;
	
	std::vector<mpz_class> result[2];
	
	result[0].resize(total);
	result[1].resize(total);
	for(int i = 0; i < total; ++i)
		result[0][i] = input[i] % p;
	int current;
	int nxt;
	for(int i = 0; i < NumOfVar; ++i)
	{
		current = i & 1;
		nxt = current ^ 1;
		int group_size = 1 << (i + 1);
		int half_group_size = 1 << i;
		for(int j = 0; j < (1 << NumOfVar) / group_size; ++j)
		{
			int base = j * group_size;
			for(int k = 0; k < group_size / 2; ++k)
			{
				result[nxt][base + k] = result[current][base + half_group_size + k];
				result[nxt][base + k + half_group_size] = (result[current][base + k] - result[current][base + half_group_size + k] + p) % p;
			}
		}
	}
	return result[nxt];
}

void test(int l){
	clock_t preinput_t = clock();
	NumOfVar = l;
	std::vector<mpz_class> test(pow(2, l));
	for(int i = 0; i < pow(2, l); i++)
		test[i] = rand();
	std::vector<mpz_class> res = pre_input(test);
	//for(int i = 0; i < 16; i++)
		//cout << "res[i] = " << res[i] << endl;
	cout << "pre_input time: " << (double)(clock() - preinput_t) / CLOCKS_PER_SEC << endl;
	return;
}	

//std::vector<Ec1> g1_pre;
void precompute_g1(){
	g1_pre.resize(P);
	g2_pre.resize(P);
	g1_pre[0] = g1;
	g2_pre[0] = g2;
	for(int i = 1; i < P; i++){
		g1_pre[i] = g1_pre[i-1] + g1_pre[i-1];
		g2_pre[i] = g2_pre[i-1] + g2_pre[i-1];
	}
	return;
}
/*
Ec1 multi_scalar(vector<Ec1>& pub, vector<Ec1>& pre, vector<mpz_class>& e){
	Ec1 temp = pubs_g1[0][0]*0;
	
	clock_t t1=clock();
	
	int length = P;

	for(int i=length-1;i>=0;i--){
		for(int j=0;j<e.size()/W;j++){
			int selector = 0;
			for(int k=W-1;k>=1;k--){
				selector+=mpz_tstbit(e[j*W+k].get_mpz_t(),i);
				selector*=2;
				
			}
			selector+=mpz_tstbit(e[j*W].get_mpz_t(),i);

			if(selector == 3) 
				temp+=pre[j];
			else if(selector == 1)
				temp+=pub[j*2];
			else if(selector == 2)
				temp+=pub[j*2+1];
			else{}	
			
		}
		if(i!=0)
			temp = temp*2;
			
	}
	return temp;
}
*/
template<class T>
T pre_exp(vector<T>& pre, mpz_class n){
	T temp = pre[0]*0;
	int length = mpz_sizeinbase(n.get_mpz_t(), 2);
	//cout << "length = " << length << endl;
	for(int i = 0; i < length; i++){
		if(mpz_tstbit(n.get_mpz_t(),i) == 1)
			temp = temp + pre[i];
	}
	return temp;
}

Ec1 g_baby_step[256 / 15 + 1][1 << 15];

void KeyGen_preprocessing(Ec1 g)
{
	printf("Preprocess start\n");
	Ec1 g_pow = g;
	mie::Vuint exponent = 1;
	for(int i = 0; i < 256 / 15 + 1; ++i)
	{
		g_baby_step[i][0] = g * 0;
		for(int j = 1; j < (1 << 15); ++j)
			g_baby_step[i][j] = g_baby_step[i][j - 1] + g_pow;
		g_pow = g_pow * (1 << 15);
	}
	printf("Preprocess end\n");
}

Ec1 g1_exp(mpz_class a)
{
	Ec1 ret = g1 * 0;
	int length = mpz_sizeinbase(a.get_mpz_t(), 2);
	for(int i = 0; i < length / 15 + 1; ++i)
	{
		mpz_class mask_mpz = (a >> (i * 15)) & ((1 << 15) - 1);
		int mask = mpz_get_si(mask_mpz.get_mpz_t());
		ret = ret + g_baby_step[i][mask];
	}
	return ret;
}

void KeyGen(int d){
	NumOfVar = d;
	clock_t KeyGen_t = clock();
	mpz_urandomm(a.get_mpz_t(), r_state, p.get_mpz_t());
	
	//const mie::Vuint temp(a.get_str().c_str());
	
	//g1a = g1 * temp;
	//g2a = g2 * temp;
	precompute_g1();
	KeyGen_preprocessing(g1);
	g1a = pre_exp(g1_pre, a);
	g2a = pre_exp(g2_pre, a);

	s.resize(NumOfVar);
	for(int i = 0; i < NumOfVar; i++)
		mpz_urandomm(s[i].get_mpz_t(), r_state, p.get_mpz_t());
	pub_g1_exp.resize((int)pow(2, d));
	pub_g1.resize((int)pow(2, d));
	pub_g2.resize((int)pow(2, d));
	pub_g1_exp[0] = 1;
	pub_g1[0] = g1;
	pub_g2[0] = g2;
	//precompute_g1();
	for(int i = 0; i < d; i++){
		for(int j = 1 << i; j < (1 << (i + 1)); j++){
			pub_g1_exp[j] = (s[i] * pub_g1_exp[j - (1 << i)]) % p;
			pub_g1[j] = g1_exp(pub_g1_exp[j]);
				//pub_g2[j] = g2 * temp1;
		}
	}

	//multi_scalar
	vector<Ec1> scalars;
	scalars.resize(multi_scalar_w);
	multi_scalar_g1.resize((1 << d) / multi_scalar_w + 1);
	for(int i = 0; i < (1 << d) / multi_scalar_w + 1; ++i)
	{
		for(int j = 0; j < multi_scalar_w; ++j)
		{
			int id = i * multi_scalar_w + j;
			if(j >= (1 << d))
			{
				scalars[j] = g1 * 0;
			}
			else
			{
				scalars[j] = pub_g1[id];
			}
		}
		for(int j = 0; j < (1 << multi_scalar_w); ++j)
		{
			multi_scalar_g1[i].value[j] = g1 * 0;
			for(int k = 0; k < multi_scalar_w; ++k)
			{
				if((j >> k) & 1)
				{
					multi_scalar_g1[i].value[j] = multi_scalar_g1[i].value[j] + scalars[k];
				}
			}
		}
	}
	for(int i = 0; i < d; ++i)
	{
		pub_g2[1 << i] = pre_exp(g2_pre, s[i]);
	}
	cout << "KeyGen time: " << (double)(clock() - KeyGen_t) / CLOCKS_PER_SEC << endl;
	
	return;
}

Ec1 multi_scalar_calc(int index, const vector<mpz_class> &scalar_pow)
{
	Ec1 ret = g1 * 0;
	int max_len = 0;
	for(int j = 0; j < multi_scalar_w; ++j)
		max_len = max(max_len, (int)mpz_sizeinbase(scalar_pow[j].get_mpz_t(), 2));
	for(int j = max_len - 1; j >= 0; --j)
	{
		int current_scalar = 0;
		for(int k = 0; k < multi_scalar_w; ++k)
			current_scalar |= (mpz_tstbit(scalar_pow[k].get_mpz_t(), j) == 1) << (k);
		ret = ret * 2 + multi_scalar_g1[index].value[current_scalar];
	}
	return ret;
}

void commit(Ec1& digest, Ec1& digesta, vector<mpz_class>& input){
	
	vector<mpz_class> coeffs = input;
	
	clock_t commit_t = clock();
	
	//compute digest pub
	
	for(int i = 0; i < coeffs.size(); i++)
		if(coeffs[i] < 0)
			coeffs[i] += p;
	//vector<Ec1> pub_pre(2 * NumOfVar + 1);
	//mpz_class ans = 0;
	vector<mpz_class> scalar_pow;
	scalar_pow.resize(multi_scalar_w);

	for(int i = 0; i < (1 << NumOfVar) / multi_scalar_w + 1; i++){
		for(int j = 0; j < multi_scalar_w; ++j)
		{
			int id = i * multi_scalar_w + j;
			if(id >= (1 << NumOfVar))
			{
				scalar_pow[j] = 0;
			}
			else
			{
				scalar_pow[j] = coeffs[id];
			}
		}
		//cout << "i = " << i << endl;
		//cout << "pub_g1[i] * temp = " << pub_g1[i] * temp << endl;
		digest = digest + multi_scalar_calc(i, scalar_pow);
	} 
	const mie::Vuint temp1(a.get_str().c_str());

	digesta = digest * temp1;
	
	cout << "commit time: " << (double)(clock() - commit_t) / CLOCKS_PER_SEC << endl;
	
	return;
	
}

bool check_commit(Ec1 digest, Ec1 digesta){
	Fp12 ea1, ea2;
	
	opt_atePairing(ea1, g2, digesta);
	opt_atePairing(ea2, g2a, digest);
	
	
	return (ea1 == ea2);
}

void prove(vector<mpz_class> r, mpz_class& ans, vector<mpz_class>& input, vector<Ec1>& witness, vector<Ec1>& witnessa){
	vector<mpz_class> coeffs = input;
	clock_t prove_t = clock();

	for(int i = 0; i < coeffs.size(); i++)
		if(coeffs[i] < 0)
			coeffs[i] += p;
	//vector<Ec1> pub_pre(2 * NumOfVar + 1);
	//cout << "ans = " << ans << endl;
	std::vector<mpz_class> ans_pre;
	ans_pre.resize((int)pow(2, NumOfVar));
	ans_pre[0] = 1;
	ans = coeffs[0];
	for(int i = 0; i < NumOfVar; i++){
		for(int j = (int)pow(2, i); j < (int)pow(2, i+1); j++){
			//cout << "j = " << j << endl;
			ans_pre[j] = (r[i] * ans_pre[j - (int)pow(2, i)]) % p;
			ans = (ans + (ans_pre[j] * coeffs[j]) % p) % p; 
		}	
	}
	//cout << "ans = " << ans << endl;
	//cout << "rightans = " << (coeffs[0] + r[0] * coeffs[1] + r[1] * coeffs[2] + r[0] * r[1] * coeffs[3]) % p << endl;

	witness.resize(NumOfVar);
	for(int i = 0; i < NumOfVar; i++)
		witness[i] = g1 * 0;
	witnessa.resize(NumOfVar);
	vector<mpz_class> witness_coeffs((int)pow(2, NumOfVar)), temp_coeffs = coeffs;

	int start_index = 0;
	
	for(int i = 0; i < NumOfVar; i++){
		for(int j = 0; j < pow(2,r.size() - i - 1); j++){
			witness_coeffs[start_index + j] = temp_coeffs[2 * j + 1] % p;
			temp_coeffs[j] = (temp_coeffs[2 * j] + temp_coeffs[2 * j + 1] * r[i]) % p;
			//cout << "i = " << i << " j = " << j << " witness = " << witness_coeffs[start_index + j] << endl;
		}
		temp_coeffs.resize(pow(2, r.size() - i - 1));
		start_index += pow(2, r.size() - i - 1);
	}

	for(int i = 0; i < witness_coeffs.size(); i++)
		if(witness_coeffs[i] < 0)
			witness_coeffs[i] += p;

	/*
	int start_index2 = 0;
	for(int i = 0; i < NumOfVar; i++){
		for(int j = 0; j < pow(2, r.size() - i - 1); j++){
			mie::Vuint temp(witness_coeffs[start_index2 + j].get_str().c_str());
			witness[i] += pub_g1[pow(2, i + 1) * j] * temp; 
			mie::Vuint temp1(coeffs[1].get_str().c_str());

			//cout << "g1*coeffs[1] = " << g1 * temp1 << endl;
			//cout << "witness[0] = " << witness[i] << endl;
		}
		const mie::Vuint temp2(a.get_str().c_str());
		witnessa[i] = witness[i] * temp2;
		start_index2 += pow(2, r.size() - i - 1);
	}
	*/
	int start_index2 = 0;
	std::vector<mpz_class> witness_pre(r.size());
	for(int i = 0; i < NumOfVar; i++){
		for(int j = 0; j < pow(2, r.size() - i - 1); j++){
			//mie::Vuint temp(witness_coeffs[start_index2 + j].get_str().c_str());
			witness_pre[i] += (pub_g1_exp[pow(2, i + 1) * j] * witness_coeffs[start_index2 + j]) % p;
			//witness[i] = pre_exp(g1_pre, witness_pre[i])
			//cout << "g1*coeffs[1] = " << g1 * temp1 << endl;
			//cout << "witness[0] = " << witness[i] << endl;
		}
		witness[i] = pre_exp(g1_pre, witness_pre[i]);
		const mie::Vuint temp2(a.get_str().c_str());
		witnessa[i] = witness[i] * temp2;
		start_index2 += pow(2, r.size() - i - 1);
	}
	cout << "prove time: " << (double)(clock() - prove_t) / CLOCKS_PER_SEC << endl;	
}

bool verify(vector<mpz_class> r, Ec1 digest, mpz_class& ans, vector<Ec1>& witness, vector<Ec1>& witnessa){
	clock_t verify_t = clock();	

	Fp12 ea1, ea2;
	
	bool flag = 1;

	cout << r.size() << endl;
	
	for(int i = 0; i < r.size(); i++){
		opt_atePairing(ea1, g2, witnessa[i]);
		opt_atePairing(ea2, g2a, witness[i]);
		if(ea1 != ea2){
			cout << "here error!" << endl;
			flag = 0;
		}
	}
	
	Fp12 ea3, ea4 = 1;

	std::vector<Fp12> temp(r.size());
	//ans = ans % p;
	//globalleft = globalleft % p + p;
	//mie::Vuint temp1(ans.get_str().c_str());
	//cout << "testans = " << ans << endl;

	//Ec1 temp2 = g1 * temp1;
	Ec1 temp2 = pre_exp(g1_pre, ans);

	opt_atePairing(ea3, g2, digest - temp2);

	for(int i = 0; i < r.size(); i++){

		mie::Vuint temp4(r[i].get_str().c_str());
		
		//Ec2 temp3 = pub_g2[pow(2, i)] - g2 * temp4;
		Ec2 temp3 = pub_g2[pow(2, i)] - pre_exp(g2_pre, r[i]);
		
		opt_atePairing(temp[i], temp3, witness[i]);
		
		ea4 *= temp[i];
	}





	if(ea3 != ea4) {
		cout << "final error" << endl;
		flag = 0;
	}

	cout << "verify time: "<<(double)(clock() - verify_t) / CLOCKS_PER_SEC << endl;
	return flag;
	
}


int main(int argc, char** argv){
	p.set_str("16798108731015832284940804142231733909759579603404752749028378864165570215949",10);
	test(atoi(argv[1]));
	seed = rand();
    gmp_randinit_default(r_state);
    gmp_randseed_ui(r_state, seed);
	
	
	//bilinear g1 g2
	bn::CurveParam cp = bn::CurveFp254BNb;
	Param::init(cp);
	const Point& pt = selectPoint(cp);
	g2 = Ec2(
		Fp2(Fp(pt.g2.aa), Fp(pt.g2.ab)),
		Fp2(Fp(pt.g2.ba), Fp(pt.g2.bb))
	);
	g1 = Ec1(pt.g1.a, pt.g1.b);
	

	int d = atoi(argv[1]);
	
	KeyGen(d);
	
	int N = (int)pow(2, d);
	vector<mpz_class> input(N);
	for(int i = 0; i < input.size(); i++){
		input[i] = rand();
		//cout << "input[i] = " << input[i] << endl;
	}
	
	
	
	Ec1 digest, digesta;
	digest = g1 * 0;
	
	commit(digest,digesta,input);
	cout << "check commit: " << check_commit(digest, digesta) << endl;
	
	vector<Ec1> proof, proofa;
	vector<mpz_class> r(d);
	mpz_class ans;
	
	for(int i = 0; i < r.size(); i++)
		mpz_urandomm(r[i].get_mpz_t(), r_state, p.get_mpz_t());
	
	prove(r, ans, input, proof, proofa);
	
	bool tf = verify(r, digest, ans, proof, proofa);
	cout<< "verify: " << tf <<endl;
	return 0;
}
