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

namespace vpd_test
{

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
	s.resize(NumOfVar + 1);
	for(int i = 0; i < NumOfVar + 1; i++)
		mpz_urandomm(s.at(i).get_mpz_t(), r_state, p.get_mpz_t());
	//s[NumOfVar] = 0;
	//+6 for s_{un}^5, s_{un}^4, s_{un}^3 and s_{vn}^5, s_{vn}^4, s_{vn}^3
	pub_g1.resize(d * 2 + 1 + 6 + 1);
	pub_g2.resize(d * 2 + 1 + 6 + 1);
	for(int i = 0; i < d; i++){
	
		pub_g1.at(2 * i + 1) = pre_exp(g1_pre, s.at(i));
		//pub_g2[2 * i + 1] = g2 * temp2;
		pub_g2.at(2 * i + 1) = pre_exp(g2_pre, s.at(i));
		mpz_class square = (s.at(i) * s.at(i)) % p;
		//mie::Vuint temp1(square.get_str().c_str());
		//pub_g1[2 * i] = g1 * temp1; 
		pub_g1.at(2 * i) = pre_exp(g1_pre, square);
	
		//pub_g2[2 * i] = g2 * temp1; 
		//pub_g2[2 * i] = pre_exp(g2_pre, square); 
		if(i == d / 2 - 1){
			//cout << "i = " << i << endl;
			mpz_class cubic = (square * s[i]) % p;
			mpz_class qudruple = (cubic * s[i]) % p;
			mpz_class quintuple = (qudruple * s[i]) % p;
			pub_g1.at(2 * d + 1) = pre_exp(g1_pre, quintuple);
			pub_g1.at(2 * d + 2) = pre_exp(g1_pre, qudruple);
			pub_g1.at(2 * d + 3) = pre_exp(g1_pre, cubic);
			//pub_g2[2 * d + 1] = pre_exp(g2_pre, quintuple);
			//pub_g2[2 * d + 2] = pre_exp(g2_pre, qudruple);
			//pub_g2[2 * d + 3] = pre_exp(g2_pre, cubic);
		}
		if(i == d - 2){
			//cout << "i = " << i << endl;
			mpz_class cubic = (square * s[i]) % p;
			mpz_class qudruple = (cubic * s[i]) % p;
			mpz_class quintuple = (qudruple * s[i]) % p;
			pub_g1.at(2 * d + 4) = pre_exp(g1_pre, quintuple);
			pub_g1.at(2 * d + 5) = pre_exp(g1_pre, qudruple);
			pub_g1.at(2 * d + 6) = pre_exp(g1_pre, cubic);
			//pub_g2[2 * d + 4] = pre_exp(g2_pre, quintuple);
			//pub_g2[2 * d + 5] = pre_exp(g2_pre, qudruple);
			//pub_g2[2 * d + 6] = pre_exp(g2_pre, cubic);
		}
	}
	pub_g1.at(2 * d) = g1;
	pub_g2.at(2 * d) = g2;
	pub_g1.at(2 * d + 7) = pre_exp(g1_pre, s[d]);
	std::cout << "pub_g1[2 * d + 7] = " << pub_g1[2 * d + 7] << std::endl;
	pub_g2.at(2 * d + 1) = pre_exp(g2_pre, s[d]);

	cout << "KeyGen time: " << (double)(clock() - KeyGen_t) / CLOCKS_PER_SEC << endl;
	
	return;
}

mpz_class commit(Ec1& digest, Ec1& digesta, vector<mpz_class>& input){
	//cout << "digest = " << digest << endl; 
	digest = g1 * 0;
	mpz_class r_f;
	mpz_urandomm(r_f.get_mpz_t(), r_state, p.get_mpz_t());
	vector<mpz_class> coeffs = input;
	coeffs.push_back(r_f);
	
	//int d = ceil(log2(input.size()));
	
	clock_t commit_t = clock();
	
	//compute digest pub
	
	for(int i = 0; i < coeffs.size(); i++)
		if(coeffs[i] < 0)
			coeffs[i] += p;
	//vector<Ec1> pub_pre(2 * NumOfVar + 1 + 6);
	mpz_class ans = 0;
	for(int i = 0; i < 2 * NumOfVar + 1 + 6 + 1; i++){
		mie::Vuint temp(coeffs.at(i).get_str().c_str());
		digest = digest + (pub_g1.at(i) * temp);
	} 

	mie::Vuint temp1(a.get_str().c_str());

	digesta = digest * temp1;
	
	//result = digest;
	
	cout << "commit time: " << (double)(clock() - commit_t) / CLOCKS_PER_SEC << endl;
	
	return r_f;
	
}

bool check_commit(Ec1 digest, Ec1 digesta){
	Fp12 ea1, ea2;
	
	opt_atePairing(ea1, g2, digesta);
	opt_atePairing(ea2, g2a, digest);
	
	
	return (ea1 == ea2);
}

void prove(vector<mpz_class> r, mpz_class& ans, vector<mpz_class>& input, vector<Ec1>& witness, vector<Ec1>& witnessa, mpz_class r_f){
	vector<mpz_class> coeffs = input;
	clock_t prove_t = clock();

	witness.resize(NumOfVar + 1);
	witnessa.resize(NumOfVar + 1);
	std::vector<mpz_class> t(NumOfVar);
	for(int i = 0; i < t.size(); i++)
		mpz_urandomm(t[i].get_mpz_t(), r_state, p.get_mpz_t());
	coeffs.push_back(r_f);
	for(int i = 0; i < NumOfVar; i++){
		if(i == NumOfVar / 2 - 1){
			mpz_class tmp1 = coeffs[2 * NumOfVar + 1];
			mie::Vuint temp1(tmp1.get_str().c_str());
			mpz_class tmp2 = (tmp1 * r[i] + coeffs[2 * NumOfVar + 2]) % p;
			mie::Vuint temp2(tmp2.get_str().c_str());
			mpz_class tmp3 = (tmp2 * r[i] + coeffs[2 * NumOfVar + 3]) % p;
			mie::Vuint temp3(tmp3.get_str().c_str());
			mpz_class tmp4 = (tmp3 * r[i] + coeffs[2 * i]) % p;
			mie::Vuint temp4(tmp4.get_str().c_str());
			mpz_class tmp5 = (tmp4 * r[i] + coeffs[2 * i + 1]) % p;
			mie::Vuint temp5(tmp5.get_str().c_str());

			witness[i] = pub_g1[2 * NumOfVar + 2] * temp1 + pub_g1[2 * NumOfVar + 3] * temp2 + pub_g1[2 * i] * temp3 + pub_g1[2 * i + 1] * temp4 + pre_exp(g1_pre, tmp5);
		}
		if(i == NumOfVar - 2){
			mpz_class tmp1 = coeffs[2 * NumOfVar + 4];
			mie::Vuint temp1(tmp1.get_str().c_str());
			mpz_class tmp2 = (tmp1 * r[i] + coeffs[2 * NumOfVar + 5]) % p;
			mie::Vuint temp2(tmp2.get_str().c_str());
			mpz_class tmp3 = (tmp2 * r[i] + coeffs[2 * NumOfVar + 6]) % p;
			mie::Vuint temp3(tmp3.get_str().c_str());
			mpz_class tmp4 = (tmp3 * r[i] + coeffs[2 * i]) % p;
			mie::Vuint temp4(tmp4.get_str().c_str());
			mpz_class tmp5 = (tmp4 * r[i] + coeffs[2 * i + 1]) % p;
			mie::Vuint temp5(tmp5.get_str().c_str());

			witness[i] = pub_g1[2 * NumOfVar + 5] * temp1 + pub_g1[2 * NumOfVar + 6] * temp2 + pub_g1[2 * i] * temp3 + pub_g1[2 * i + 1] * temp4 + pre_exp(g1_pre, tmp5);
		}
		if(i != NumOfVar / 2 - 1 && i != NumOfVar - 2){
			mie::Vuint temp1(coeffs[2 * i].get_str().c_str());
			mpz_class tmp = (coeffs[2 * i] * r[i]) % p;
			if(tmp < 0) tmp += p;
			
			witness[i] = pub_g1[2 * i + 1] * temp1 + pre_exp(g1_pre, tmp) + pre_exp(g1_pre, coeffs[2 * i + 1]);
		}
		mie::Vuint tempti(t[i].get_str().c_str());
		witness[i] += pub_g1[2 * NumOfVar + 7] * tempti;
		mie::Vuint tempa(a.get_str().c_str());
		witnessa[i] = witness[i] * tempa;
	}
	mpz_class tmp = coeffs[2 * NumOfVar + 7];
	for(int i = 0; i < NumOfVar; i++)
		tmp += (t[i] * r[i]) % p;
	//mie::Vuint temp(tmp.get_str().c_str());
	witness[NumOfVar] = pre_exp(g1_pre, tmp);
	clock_t zkt = clock();
	for(int i = 0; i < NumOfVar; i++){
		mie::Vuint tempti(t[i].get_str().c_str());
		witness[NumOfVar] -= pub_g1[2 * i + 1] * tempti;
	}
	mie::Vuint tempa(a.get_str().c_str());
	witnessa[NumOfVar] = witness[NumOfVar] * tempa;
	cout << "zkt time: " << (double)(clock() - zkt) / CLOCKS_PER_SEC << endl;	

	for(int i = 0; i < coeffs.size(); i++)
		if(coeffs[i] < 0)
			coeffs[i] += p;
	//vector<Ec1> pub_pre(2 * NumOfVar + 1);
	//cout << "ans = " << ans << endl;
	ans = 0;
	for(int i = 0; i < NumOfVar; i++){
		if(i == NumOfVar / 2 - 1){
			mpz_class tmp = coeffs[2 * NumOfVar + 1];
			tmp = (tmp * r[i] + coeffs[2 * NumOfVar + 2]) % p;
			tmp = (tmp * r[i] + coeffs[2 * NumOfVar + 3]) % p;
			tmp = (tmp * r[i] + coeffs[2 * i]) % p;
			tmp = (tmp * r[i] + coeffs[2 * i + 1]) % p;
			tmp = (tmp * r[i]) % p;
			ans = (ans + tmp) % p;
		}
		if(i == NumOfVar - 2){
			mpz_class tmp = coeffs[2 * NumOfVar + 4];
			tmp = (tmp * r[i] + coeffs[2 * NumOfVar + 5]) % p;
			tmp = (tmp * r[i] + coeffs[2 * NumOfVar + 6]) % p;
			tmp = (tmp * r[i] + coeffs[2 * i]) % p;
			tmp = (tmp * r[i] + coeffs[2 * i + 1]) % p;
			tmp = (tmp * r[i]) % p;
			ans = (ans + tmp) % p;
		}
		if(i != NumOfVar / 2 - 1 && i != NumOfVar - 2){
			ans = (ans + coeffs[2 * i] * r[i] * r[i]) % p;
			//cout << "ans = " << ans << endl;
			ans = (ans + coeffs[2 * i + 1] * r[i]) % p;
		}
		//cout << "ans = " << ans << endl;
	}
	ans = (ans + coeffs[2 * NumOfVar]) % p;
	cout << "prove time: " << (double)(clock() - prove_t) / CLOCKS_PER_SEC << endl;	
}

bool verify(vector<mpz_class> r, Ec1 digest, mpz_class& ans, vector<Ec1>& witness, vector<Ec1>& witnessa){
	clock_t verify_t = clock();	

	Fp12 ea1, ea2;
	
	bool flag = 1;
	
	for(int i = 0; i < r.size() + 1; i++){
		opt_atePairing(ea1, g2, witnessa[i]);
		opt_atePairing(ea2, g2a, witness[i]);
		if(ea1 != ea2){
			cout << "here error!" << endl;
			flag = 0;
		}
	}
	
	Fp12 ea3, ea4 = 1;

	std::vector<Fp12> temp(r.size() + 1);
	ans = ans % p;
	Ec1 temp2 = pre_exp(g1_pre, ans);

	opt_atePairing(ea3, g2, digest - temp2);

	for(int i = 0; i < r.size() + 1; i++){
		//cout << "i = " << i << endl;
		if(i < r.size()){
			Ec2 temp3 = pub_g2[2 * i + 1] - pre_exp(g2_pre, r[i]);
			opt_atePairing(temp[i], temp3, witness[i]);
		}
		if(i == r.size()){
			//cout << "i = " << i << endl;
			opt_atePairing(temp[i], pub_g2[2 * i + 1], witness[i]);
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

void environment_init()
{
	cout << "vpd_test initialized" << endl;
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
}

}
