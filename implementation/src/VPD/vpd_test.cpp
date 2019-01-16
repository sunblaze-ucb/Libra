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
	s.resize(NumOfVar);
	for(int i = 0; i < NumOfVar; i++)
		mpz_urandomm(s[i].get_mpz_t(), r_state, p.get_mpz_t());
	pub_g1.resize(d * 2 + 1);
	pub_g2.resize(d * 2 + 1);
	for(int i = 0; i < d; i++){
		mpz_class square = (s[i] * s[i]) % p;
		//mie::Vuint temp1(square.get_str().c_str());
		//pub_g1[2 * i] = g1 * temp1; 
		pub_g1[2 * i] = pre_exp(g1_pre, square);
	
		//pub_g2[2 * i] = g2 * temp1; 
		pub_g2[2 * i] = pre_exp(g2_pre, square);  
	}
	for(int i = 0; i < d; i++){
		//mie::Vuint temp2(s[i].get_str().c_str());
		//pub_g1[2 * i + 1] = g1 * temp2; 
		pub_g1[2 * i + 1] = pre_exp(g1_pre, s[i]);
		//pub_g2[2 * i + 1] = g2 * temp2;
		pub_g2[2 * i + 1] = pre_exp(g2_pre, s[i]);
	}
	pub_g1[2 * d] = g1;
	pub_g2[2 * d] = g2;
	cout << "KeyGen time: " << (double)(clock() - KeyGen_t) / CLOCKS_PER_SEC << endl;
	
	return;
}

void commit(Ec1& digest, Ec1& digesta, vector<mpz_class>& input){
	//cout << "digest = " << digest << endl; 
	vector<mpz_class> coeffs = input;
	
	//int d = ceil(log2(input.size()));
	
	clock_t commit_t = clock();
	
	//compute digest pub
	
	for(int i = 0; i < coeffs.size(); i++)
		if(coeffs[i] < 0)
			coeffs[i] += p;
	vector<Ec1> pub_pre(2 * NumOfVar + 1);
	mpz_class ans = 0;
	for(int i = 0; i < 2 * NumOfVar + 1; i++){
		mie::Vuint temp(coeffs[i].get_str().c_str());
		digest = digest + (pub_g1[i] * temp);
	} 

	mie::Vuint temp1(a.get_str().c_str());

	digesta = digest * temp1;
	
	result = digest;
	
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

	witness.resize(NumOfVar);
	witnessa.resize(NumOfVar);
	for(int i = 0; i < NumOfVar; i++){
		mie::Vuint temp1(coeffs[2 * i].get_str().c_str());
		mpz_class tmp = (coeffs[2 * i] * r[i]) % p;
		if(tmp < 0) tmp += p;
		//mie::Vuint temp2(tmp.get_str().c_str());
		//mie::Vuint temp3(coeffs[2 * i + 1].get_str().c_str());
		//witness[i] = pub_g1[2 * i + 1] * temp1 + g1 * temp2 + g1 * temp3;
		witness[i] = pub_g1[2 * i + 1] * temp1 + pre_exp(g1_pre, tmp) + pre_exp(g1_pre, coeffs[2 * i + 1]);
		mie::Vuint temp4(a.get_str().c_str());
		witnessa[i] = witness[i] * temp4;
	}
	for(int i = 0; i < coeffs.size(); i++)
		if(coeffs[i] < 0)
			coeffs[i] += p;
	//vector<Ec1> pub_pre(2 * NumOfVar + 1);
	//cout << "ans = " << ans << endl;
	for(int i = 0; i < NumOfVar; i++){
		ans = ans + coeffs[2 * i] * r[i] * r[i];
		//cout << "ans = " << ans << endl;
		ans = ans + coeffs[2 * i + 1] * r[i];
		//cout << "ans = " << ans << endl;
	}
	ans = ans + coeffs[2 * NumOfVar];
	
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

	std::vector<Fp12> temp(r.size());
	ans = ans % p;
	//globalleft = globalleft % p + p;
	mie::Vuint temp1(ans.get_str().c_str());
	//cout << "testans = " << ans << endl;

	//Ec1 temp2 = g1 * temp1;
	Ec1 temp2 = pre_exp(g1_pre, ans);

	opt_atePairing(ea3, g2, digest - temp2);

	for(int i = 0; i < r.size(); i++){
		//mpz_class sr = (s[i] - r[i]) % p;
		//if(sr < 0) sr += p; 
		mie::Vuint temp4(r[i].get_str().c_str());
		//cout << "sr = " << sr << endl;
		//cout << "temp4 = " << temp4 << endl;
		//Ec2 temp3 = pub_g2[2 * i + 1] - g2 * temp4;
		Ec2 temp3 = pub_g2[2 * i + 1] - pre_exp(g2_pre, r[i]);
		opt_atePairing(temp[i], temp3, witness[i]);
		//cout << "globalleft = " << globalleft << endl;
		//cout << "globalright = " << globalright * sr << endl;
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
	

	int d = atoi(argv[1]);
	
	KeyGen(d);
	
	int N = NumOfVar * 2 + 1;
	vector<mpz_class> input(N);
	for(int i = 0; i < input.size(); i++){
		input[i] = rand();
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



