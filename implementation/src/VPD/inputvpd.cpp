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

#include "mcl/bn256.hpp"

using namespace std;
using namespace mcl::bn256;

#define P 512
int multi_scalar_w;

unsigned long int seed;
int NumOfVar;
gmp_randstate_t r_state;  
mpz_class p, a;

G1 g1, g1a, result;
G2 g2, g2a;

vector<G1> pub_g1, g1_pre;
vector<G2> pub_g2, g2_pre;

class multi_scalar_state
{
public:
<<<<<<< HEAD
<<<<<<< HEAD
	G1 value[1 << multi_scalar_w];
=======
	//Ec1 value[1 << multi_scalar_w];
	Ec1 value[1 << 10];
>>>>>>> parent of 8c81323... minor
=======
	//Ec1 value[1 << multi_scalar_w];
	Ec1 value[1 << 10];
>>>>>>> parent of 8c81323... minor
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
	//cout << "pre_input time: " << (double)(clock() - preinput_t) / CLOCKS_PER_SEC << endl;
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

G1 g_baby_step[256 / 15 + 1][1 << 15];

void KeyGen_preprocessing(G1 g)
{
	//printf("Preprocess start\n");
	G1 g_pow = g;
	for(int i = 0; i < 256 / 15 + 1; ++i)
	{
		g_baby_step[i][0] = g * 0;
		for(int j = 1; j < (1 << 15); ++j)
			g_baby_step[i][j] = g_baby_step[i][j - 1] + g_pow;
		g_pow = g_pow * (1 << 15);
	}
	//printf("Preprocess end\n");
}

G1 g1_exp(mpz_class a)
{
	G1 ret = g1 * 0;
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

	s.resize(NumOfVar + 1);
	for(int i = 0; i < NumOfVar + 1; i++)
		mpz_urandomm(s[i].get_mpz_t(), r_state, p.get_mpz_t());
	pub_g1_exp.resize((int)pow(2, d));
	pub_g1.resize((int)pow(2, d) + 1);
	pub_g2.resize(d + 1);
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
	pub_g1[1 << d] = pre_exp(g1_pre, s[d]);
	//multi_scalar
	vector<G1> scalars;
	scalars.resize(multi_scalar_w);
	multi_scalar_g1.resize((1 << d) / multi_scalar_w + 1);
	for(int i = 0; i < (1 << d) / multi_scalar_w + 1; ++i)
	{
		for(int j = 0; j < multi_scalar_w; ++j)
		{
			int id = i * multi_scalar_w + j;
			if(id >= (1 << d))
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
	for(int i = 0; i < d + 1; ++i)
	{
		pub_g2[i] = pre_exp(g2_pre, s[i]);
	}
	//cout << "KeyGen time: " << (double)(clock() - KeyGen_t) / CLOCKS_PER_SEC << endl;
	
	return;
}

G1 multi_scalar_calc(int index, const vector<mpz_class> &scalar_pow)
{
	G1 ret = g1 * 0;
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

void commit(G1& digest, G1& digesta, vector<mpz_class>& input){
	
	vector<mpz_class> coeffs = input;
	
	clock_t commit_t = clock();
	
	//compute digest pub
	
	for(int i = 0; i < coeffs.size(); i++)
		if(coeffs[i] < 0)
			coeffs[i] += p;
	//vector<G1> pub_pre(2 * NumOfVar + 1);
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

	digest += pub_g1[1 << NumOfVar] * input[1 << NumOfVar];

	digesta = digest * a;

	//cout << "commit time: " << (double)(clock() - commit_t) / CLOCKS_PER_SEC << endl;
	
	return;
	
}

bool check_commit(G1 digest, G1 digesta){
	Fp12 ea1, ea2;
	
	pairing(ea1, digesta, g2);
	pairing(ea2, digest, g2a);
	
	
	return (ea1 == ea2);
}

void prove(vector<mpz_class> r, mpz_class& ans, vector<mpz_class>& input, vector<G1>& witness, vector<G1>& witnessa){
	vector<mpz_class> coeffs = input;
	clock_t prove_t = clock();

	std::vector<mpz_class> t(NumOfVar);
	for(int i = 0; i < t.size(); i++)
		mpz_urandomm(t[i].get_mpz_t(), r_state, p.get_mpz_t());

	for(int i = 0; i < coeffs.size(); i++)
		if(coeffs[i] < 0)
			coeffs[i] += p;
	//vector<G1> pub_pre(2 * NumOfVar + 1);
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
	witness.resize(NumOfVar + 1);
	for(int i = 0; i < NumOfVar; i++)
		witness[i] = g1 * 0;
	witnessa.resize(NumOfVar + 1);
	vector<mpz_class> witness_coeffs((int)pow(2, NumOfVar)), temp_coeffs = coeffs;

	int start_index = 0;
	//cout << "r.size() = " << r.size() << endl;
	for(int i = 0; i < NumOfVar; i++){
		for(int j = 0; j < pow(2, r.size() - i - 1); j++){
			witness_coeffs[start_index + j] = temp_coeffs[pow(2, r.size() - i - 1) + j] % p;
			temp_coeffs[j] = (temp_coeffs[j] + (temp_coeffs[pow(2, r.size() - i - 1) + j] * r[NumOfVar - 1 - i]) % p) % p;
			//cout << "p = " << temp_coeffs[j] << endl;
			//cout << "i = " << i << " j = " << j << "witness_coeffs = " << witness_coeffs[j] << endl;

		}
		//temp_coeffs.resize(pow(2, r.size() - i - 1));
		start_index += pow(2, r.size() - i - 1);
	}
	for(int i = 0; i < witness_coeffs.size(); i++)
		if(witness_coeffs[i] < 0)
			witness_coeffs[i] += p;

	//std::cout << "hello world " << std::endl;
	std::vector<mpz_class> scalar_pow(multi_scalar_w);
	int temp = 0;
	for(int k = NumOfVar - 1; k >= 0; k--){
		int temp1 = 1 << k;
		//temp += temp1;
		for(int i = 0; i < temp1 / multi_scalar_w + 1; i++){
			for(int j = 0; j < multi_scalar_w; ++j)
			{
				int id = i * multi_scalar_w + j;
				if(id >= (1 << k))
				{
					scalar_pow[j] = 0;
				}
				else
				{
					//cout << "k = " << k << " i = " << i << " j = " << j << endl;
					scalar_pow[j] = witness_coeffs[id + temp];
					//cout << "witness_coeffs[id + temp] = " << witness_coeffs[id + temp] << endl;
				}
			}
			witness[k] = witness[k] + multi_scalar_calc(i, scalar_pow);
			//std::cout << "k = " << k << " i = " << i << std::endl;
		}
		//cout << "witness[k] = " << witness[k] << endl;
		witness[k] += pub_g1[1 << NumOfVar] * t[k];
		//mie::Vuint tempa(a.get_str().c_str());
		witnessa[k] = witness[k] * a;
		temp += temp1;
		//std::cout << "hello world " << std::endl;
	}	

	mpz_class tmp = coeffs[1 << NumOfVar];
	for(int i = 0; i < NumOfVar; i++)
		tmp += (t[i] * r[i]) % p;
	//mie::Vuint temp(tmp.get_str().c_str());
	witness[NumOfVar] = pre_exp(g1_pre, tmp);
	clock_t zkt = clock();
	for(int i = 0; i < NumOfVar; i++){
		witness[NumOfVar] -= pub_g1[1 << i] * t[i];
	}
	//mie::Vuint tempa(a.get_str().c_str());
	witnessa[NumOfVar] = witness[NumOfVar] * a;
	cout << "zkt time: " << (double)(clock() - zkt) / CLOCKS_PER_SEC << endl;

	cout << "prove time: " << (double)(clock() - prove_t) / CLOCKS_PER_SEC << endl;	
}

bool verify(vector<mpz_class> r, G1 digest, mpz_class& ans, vector<G1>& witness, vector<G1>& witnessa){
	clock_t verify_t = clock();	

	Fp12 ea1, ea2;
	
	bool flag = 1;

	//cout << r.size() << endl;
	
	for(int i = 0; i < r.size(); i++){
		pairing(ea1, witnessa[i], g2);
		pairing(ea2, witness[i], g2a);
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

	//G1 temp2 = g1 * temp1;
	G1 temp2 = pre_exp(g1_pre, ans);

	pairing(ea3, digest - temp2, g2);

	for(int i = 0; i < r.size() + 1; i++){
		if(i == r.size()){
			pairing(temp[i], witness[i], pub_g2[i]);
		}
		if(i < r.size()){
			G2 temp3 = pub_g2[i] - pre_exp(g2_pre, r[i]);
		
			pairing(temp[i], witness[i], temp3);
		}
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
	
	const char *aa = "12723517038133731887338407189719511622662176727675373276651903807414909099441";
	const char *ab = "4168783608814932154536427934509895782246573715297911553964171371032945126671";
	const char *ba = "13891744915211034074451795021214165905772212241412891944830863846330766296736";
	const char *bb = "7937318970632701341203597196594272556916396164729705624521405069090520231616";
	
	//bilinear g1 g2
	mcl::bn256::CurveParam cp = mcl::bn256::CurveFp254BNb;
	//Param::init(cp);
	initPairing();
	g2 = G2(Fp2(aa, ab), Fp2(ba, bb));
	g1 = G1(-1, 1);
	

	int d = atoi(argv[1]);
	for(multi_scalar_w = 2; multi_scalar_w <= 2; multi_scalar_w++){
		cout << "multi_scalar_w = " << multi_scalar_w << endl;
	
<<<<<<< HEAD
<<<<<<< HEAD
	//cout << "p = " << p << endl;

	G1 digest, digesta;
	digest = g1 * 0;
=======
=======
>>>>>>> parent of 8c81323... minor
		KeyGen(d);
		int N = (int)pow(2, d);
		vector<mpz_class> input(N + 1);
		for(int i = 0; i < input.size(); i++){
			input[i] = rand();
		}
<<<<<<< HEAD
>>>>>>> parent of 8c81323... minor
=======
>>>>>>> parent of 8c81323... minor
		
		//cout << "p = " << p << endl;

		Ec1 digest, digesta;
		digest = g1 * 0;
			
		commit(digest,digesta,input);
		//cout << "check commit: " << check_commit(digest, digesta) << endl;
			
		vector<Ec1> proof, proofa;
		vector<mpz_class> r(d);
		mpz_class ans;
		
<<<<<<< HEAD
<<<<<<< HEAD
	vector<G1> proof, proofa;
	vector<mpz_class> r(d);
	mpz_class ans;
	
	for(int i = 0; i < r.size(); i++)
		mpz_urandomm(r[i].get_mpz_t(), r_state, p.get_mpz_t());
	
	prove(r, ans, input, proof, proofa);
	
	bool tf = verify(r, digest, ans, proof, proofa);
	cout<< "verify: " << tf <<endl;
=======
=======
>>>>>>> parent of 8c81323... minor
		for(int i = 0; i < r.size(); i++)
			mpz_urandomm(r[i].get_mpz_t(), r_state, p.get_mpz_t());
		
		prove(r, ans, input, proof, proofa);
		
		bool tf = verify(r, digest, ans, proof, proofa);
		cout<< "verify: " << tf <<endl;
	}
<<<<<<< HEAD
>>>>>>> parent of 8c81323... minor
=======
>>>>>>> parent of 8c81323... minor

	return 0;
}
