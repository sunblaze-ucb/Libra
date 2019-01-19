#include <cstdlib>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <gmp.h>
#include <gmpxx.h>
#include <math.h>
#include <string>

#include "test_point.hpp"
#include "bn.h"

using namespace std;
using namespace bn;

#define P 512

unsigned long int seed;
int NumOfVar;
gmp_randstate_t r_state;  
mpz_class p, a;
Ec1 g1, g1a, result;
Ec2 g2, g2a;

vector<Ec1> pub_g1, g1_pre;
vector<Ec2> pub_g2, g2_pre;

//mpz_class globalright, globalleft;

std::vector<mpz_class> s;

void precompute_g1(mpz_class a){
	g1_pre.resize(P);
	g2_pre.resize(P);
	g1_pre[0] = g1;
	g2_pre[0] = g2;	
	for(int i = 1; i < P; i++){
		g1_pre[i] = g1_pre[i-1] + g1_pre[i-1];
		g2_pre[i] = g2_pre[i-1] + g2_pre[i-1];
	}
}

template< class T >
T pre_exp(vector<T>& pre, mpz_class n){
	T temp = pre[0] * 0;
	int length = mpz_sizeinbase(n.get_mpz_t(), 2);
	//cout << "length = " << length << endl;
	for(int i = 0; i < length; i++){
		//cout << "hello world i = " << i << endl;
		if(mpz_tstbit(n.get_mpz_t(), i) == 1){
			//cout << "i = " << i << endl; 
			temp = temp + pre[i];
			//cout << "i = " << i << endl;		
		}
	}
	//cout << "get out" << endl;
	return temp;
}

//F = a0 + a1g1 + a2g1^2 + a3c + a4c^2 + a5g1c;
//d always equals to 2.
void KeyGen(int d){
	NumOfVar = d;
	clock_t KeyGen_t = clock();
	mpz_urandomm(a.get_mpz_t(), r_state, p.get_mpz_t());
	precompute_g1(a);
	mie::Vuint temp(a.get_str().c_str());
	//cout << "a = " << a << " " << "tmep = " << temp << endl; 
	g1a = g1 * temp;
	g2a = g2 * temp;
	//cout << "g1 = " << g1 << endl;
	//cout << "g1 = " << g1 * 2 - g1 << endl;
	//vector<mpz_class> s(NumOfVar);
	s.resize(NumOfVar + 1);
	for(int i = 0; i < NumOfVar + 1; i++)
		mpz_urandomm(s[i].get_mpz_t(), r_state, p.get_mpz_t());
	//+6 for s_{un}^5, s_{un}^4, s_{un}^3 and s_{vn}^5, s_{vn}^4, s_{vn}^3
	mpz_class square0 = (s[0] * s[0]) % p;
	mpz_class square1 = (s[1] * s[1]) % p;

	pub_g1.resize(6 + 1);
	//pub_g2.resize(9);
	pub_g1[0] = g1;
	//pub_g2[0] = g2;
	pub_g1[1] = pre_exp(g1_pre, s[0]);
	//pub_g2[1] = pre_exp(g2_pre, s[0]);
	pub_g1[2] = pre_exp(g1_pre, square0);
	//pub_g2[2] = pre_exp(g2_pre, square0);
	pub_g1[3] = pre_exp(g1_pre, s[1]);
	//pub_g2[3] = pre_exp(g2_pre, s[1]);
	pub_g1[4] = pre_exp(g1_pre, square1);
	//pub_g2[4] = pre_exp(g2_pre, square1);
	pub_g1[5] = pre_exp(g1_pre, (s[0] * s[1]) % p);
	//cout << "hello world" << endl;
	//pub_g2[5] = pre_exp(g2_pre, (s[0] * s[1]) % p);
	pub_g1[6] = pre_exp(g1_pre, s[2]);

	pub_g2.resize(2 + 1);

	pub_g2[0] = pre_exp(g2_pre, s[0]);
	pub_g2[1] = pre_exp(g2_pre, s[1]);
	pub_g2[2] = pre_exp(g2_pre, s[2]);

	cout << "KeyGen time: " << (double)(clock() - KeyGen_t) / CLOCKS_PER_SEC << endl;
	
	return;
}
//F = a0 + a1g1 + a2g1^2 + a3c + a4c^2 + a5g1c + a6g1^2c + a7g1c^2 + a8g1^2c^2;
void commit(Ec1& digest, Ec1& digesta, vector<mpz_class>& input){
	//cout << "digest = " << digest << endl; 
	vector<mpz_class> coeffs = input;
	
	//int d = ceil(log2(input.size()));
	
	clock_t commit_t = clock();
	
	//compute digest pub
	
	for(int i = 0; i < coeffs.size(); i++)
		if(coeffs[i] < 0)
			coeffs[i] += p;

	mpz_class ans = 0;
	for(int i = 0; i < coeffs.size(); i++){
		//cout << "i = " << i << endl;
		mie::Vuint temp(coeffs[i].get_str().c_str());
		digest = digest + (pub_g1[i] * temp);
	} 

	mie::Vuint temp1(a.get_str().c_str());

	digesta = digest * temp1;
	
	//result = digest;
	
	cout << "commit time: " << (double)(clock() - commit_t) / CLOCKS_PER_SEC << endl;
	
	return;
	
}

bool check_commit(Ec1 digest, Ec1 digesta){
	Fp12 ea1, ea2;
	
	opt_atePairing(ea1, g2, digesta);
	opt_atePairing(ea2, g2a, digest);
	
	
	return (ea1 == ea2);
}

//F = a0 + a1g1 + a2g1^2 + a3c + a4c^2 + a5g1c;
void prove(vector<mpz_class> r, mpz_class& ans, vector<mpz_class>& input, vector<Ec1>& witness, vector<Ec1>& witnessa){
	vector<mpz_class> coeffs = input;
	clock_t prove_t = clock();

	witness.resize(NumOfVar + 1);
	witnessa.resize(NumOfVar + 1);

	std::vector<mpz_class> t(NumOfVar);
	for(int i = 0; i < t.size(); i++)
		mpz_urandomm(t[i].get_mpz_t(), r_state, p.get_mpz_t());
	
	for(int i = 0; i < coeffs.size(); i++)
		if(coeffs[i] < 0)
			coeffs[i] += p;
	mie::Vuint tempa2(coeffs[2].get_str().c_str());
	mie::Vuint tempa5(coeffs[5].get_str().c_str());

	witness[0] = pre_exp(g1_pre, (coeffs[2] * r[0] + coeffs[1]) % p) + pub_g1[1] * tempa2 + pub_g1[3] * tempa5;
	mie::Vuint tempt0(t[0].get_str().c_str());
	witness[0] = witness[0] + pub_g1[6] * tempt0;

	mie::Vuint tempa4(coeffs[4].get_str().c_str());

	witness[1] = pre_exp(g1_pre, (coeffs[5] * r[0] + coeffs[4] * r[1] + coeffs[3]) % p) + pub_g1[3] * tempa4;
	mie::Vuint tempt1(t[1].get_str().c_str());
	witness[1] = witness[1] + pub_g1[6] * tempt1;

	mpz_class temp = (coeffs[6] + (r[0] * t[0]) % p + (r[1] * t[1]) % p) % p;
	witness[2] = pre_exp(g1_pre, temp) - pub_g1[1] * tempt0 - pub_g1[3] * tempt1; 

	for(int i = 0; i < NumOfVar + 1; i++){
		mie::Vuint tempa(a.get_str().c_str());
		witnessa[i] = witness[i] * tempa;
	}
	mpz_class square0 = (r[0] * r[0]) % p;
	mpz_class square1 = (r[1] * r[1]) % p;
	mpz_class cross = (r[0] * r[1]) % p;
	ans = coeffs[0] + ((coeffs[1] * r[0]) % p) + ((coeffs[2] * square0) % p) + ((coeffs[3] * r[1]) % p) + ((coeffs[4] * square1) % p) + ((coeffs[5] * cross) % p);
	ans = ans % p;
	if(ans < 0) ans += p;
	cout << "prove time: " << (double)(clock() - prove_t) / CLOCKS_PER_SEC << endl;	
}

bool verify(vector<mpz_class> r, Ec1 digest, mpz_class& ans, vector<Ec1>& witness, vector<Ec1>& witnessa){
	clock_t verify_t = clock();	

	Fp12 ea1, ea2;
	
	bool flag = 1;
	
	for(int i = 0; i < r.size(); i++){
		opt_atePairing(ea1, g2, witnessa[i]);
		opt_atePairing(ea2, g2a, witness[i]);
		if(ea1 != ea2){
			cout << "here error!" << endl;
			flag = 0;
		}
	}
	
	Fp12 ea3, ea4 = 1;

	std::vector<Fp12> temp(r.size() + 1);
	//ans = ans % p;
	Ec1 temp2 = pre_exp(g1_pre, ans);

	opt_atePairing(ea3, g2, digest - temp2);

	for(int i = 0; i < r.size() + 1; i++){
		//cout << "i = " << i << endl;
		if(i < r.size()){
			Ec2 temp3 = pub_g2[i] - pre_exp(g2_pre, r[i]);	
			opt_atePairing(temp[i], temp3, witness[i]);
		}
		if(i == r.size()){
			//cout << "i = " << i << endl;
			opt_atePairing(temp[i], pub_g2[i], witness[i]);
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
	seed = rand();
    gmp_randinit_default(r_state);
    gmp_randseed_ui(r_state, seed);
	p.set_str("16798108731015832284940804142231733909759579603404752749028378864165570215949",10);
	
	
	//bilinear g1 g2
	bn::CurveParam cp = bn::CurveFp254BNb;
	Param::init(cp);
	const Point& pt = selectPoint(cp);
	g2 = Ec2(
		Fp2(Fp(pt.g2.aa), Fp(pt.g2.ab)),
		Fp2(Fp(pt.g2.ba), Fp(pt.g2.bb))
	);
	g1 = Ec1(pt.g1.a, pt.g1.b);

	int d = 2;
	
//	cout << "d = " << d << endl;

	KeyGen(d);
	
	int N = 2;
	vector<mpz_class> input(6 + 1);
	for(int i = 0; i < input.size(); i++){
		input[i] = rand();
	}
	
	
	
	Ec1 digest, digesta;
	digest = g1 * 0;
	
	commit(digest, digesta, input);
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


