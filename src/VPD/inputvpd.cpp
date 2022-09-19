#include "VPD/input_vpd.h"

using namespace std;
using namespace bn;

namespace input_vpd
{
#define P 512
const int multi_scalar_w = 8;

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
	Ec1 value[(1 << multi_scalar_w) - 3];
};

vector<multi_scalar_state> multi_scalar_g1;

std::vector<mpz_class> pub_g1_exp;


vector<mpz_class> s;

void pre_input(std::vector<mpz_class>& input){
	int total = (1 << NumOfVar) + 1;
	for(int i = 0; i < (1 << (NumOfVar - 1)); ++i)
		swap(input[i], input[(1 << NumOfVar) - 1 - i]);
	for(int i = 0; i < NumOfVar; ++i)
	{
		int group_size = 1 << (i + 1);
		int half_group_size = 1 << i;
		for(int j = 0; j < (1 << NumOfVar) / group_size; ++j)
		{
			int base = j * group_size;
			for(int k = 0; k < group_size / 2; ++k)
			{
				auto tmp = input[base + k];
				input[base + k] = input[base + half_group_size + k];
				input[base + k + half_group_size] = (tmp - input[base + half_group_size + k] + p) % p;
			}
		}
	}
}

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

template<class T>
T pre_exp(vector<T>& pre, mpz_class n){
	T temp = pre[0]*0;
	int length = mpz_sizeinbase(n.get_mpz_t(), 2);
	for(int i = 0; i < length; i++){
		if(mpz_tstbit(n.get_mpz_t(),i) == 1)
			temp = temp + pre[i];
	}
	return temp;
}

Ec1 g_baby_step[256 / 15 + 1][1 << 15];

void KeyGen_preprocessing(Ec1 g)
{
	Ec1 g_pow = g;
	for(int i = 0; i < 256 / 15 + 1; ++i)
	{
		g_baby_step[i][0] = g * 0;
		for(int j = 1; j < (1 << 15); ++j)
			g_baby_step[i][j] = g_baby_step[i][j - 1] + g_pow;
		g_pow = g_pow * (1 << 15);
	}
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
	for(int i = 0; i < d; i++){
		for(int j = 1 << i; j < (1 << (i + 1)); j++){
			pub_g1_exp[j] = (s[i] * pub_g1_exp[j - (1 << i)]) % p;
			pub_g1[j] = g1_exp(pub_g1_exp[j]);
		}
	}
	pub_g1[1 << d] = pre_exp(g1_pre, s[d]);
	//multi_scalar
	//assert(multi_scalar_w == 2); //to avoid some error
	vector<Ec1> scalars;
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
		for(int j = 3; j < (1 << multi_scalar_w); ++j)
		{
			multi_scalar_g1[i].value[j - 3] = g1 * 0;
			for(int k = 0; k < multi_scalar_w; ++k)
			{
				if((j >> k) & 1)
					multi_scalar_g1[i].value[j - 3] += scalars[k];
			}
		}
	}
	for(int i = 0; i < d + 1; ++i)
	{
		pub_g2[i] = pre_exp(g2_pre, s[i]);
	}
	cout << "Input VPD KeyGen time: " << (double)(clock() - KeyGen_t) / CLOCKS_PER_SEC << endl;
	
	return;
}

Ec1 multi_scalar_calc(int index, int pub_g1_length, const vector<mpz_class> &scalar_pow)
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
		ret = ret * 2;
		switch(current_scalar)
		{
			case 0:
				break;
			case 1:
				if(index * multi_scalar_w < pub_g1_length)
					ret = ret + pub_g1[index * multi_scalar_w];
				break;
			case 2:
				if(index * multi_scalar_w + 1 < pub_g1_length)
					ret = ret + pub_g1[index * multi_scalar_w + 1];
				break;
			default:
				ret = ret + multi_scalar_g1[index].value[current_scalar - 3];
		}
	}
	return ret;
}

std::pair<mpz_class, mpz_class> commit(Ec1& digest, Ec1& digesta, Ec1& digest2, Ec1& digest2a, vector<mpz_class>& input, vector<mpz_class>& input2){

	mpz_class r_f;
	digest = g1 * 0;
	mpz_urandomm(r_f.get_mpz_t(), r_state, p.get_mpz_t());
	pre_input(input);
	vector<mpz_class>& coeffs = input;
	assert(coeffs.size() >= 1);
	coeffs[(int)coeffs.size() - 1] = r_f;
	
	clock_t commit_t = clock();
	
	//compute digest pub
	
	for(int i = 0; i < coeffs.size(); i++)
		if(coeffs[i] < 0)
			coeffs[i] += p;
	
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
		
		digest = digest + multi_scalar_calc(i, (1 << NumOfVar), scalar_pow);
	} 
	mie::Vuint temp(coeffs[1 << NumOfVar].get_str().c_str());
	digest += pub_g1[1 << NumOfVar] * temp;

	mpz_class r_f2;
	digest2 = g1 * 0;
	mpz_urandomm(r_f2.get_mpz_t(), r_state, p.get_mpz_t());

	vector<mpz_class> coeffs2 = input2;
	coeffs2.push_back(r_f2);
	for(int i = 0; i < coeffs.size(); i++)
		if(coeffs[i] < 0)
			coeffs[i] += p;
	mie::Vuint temp0(coeffs2[0].get_str().c_str());
	mie::Vuint temp1(coeffs2[1].get_str().c_str());
	mie::Vuint temprf2(coeffs2[2].get_str().c_str());
	digest2 = pub_g1[0] * temp0 + pub_g1[1 << (NumOfVar - 1)] * temp1 + pub_g1[1 << NumOfVar] * temprf2; 

	const mie::Vuint tempa(a.get_str().c_str());

	digesta = digest * tempa;
	digest2a = digest2 * tempa;

	cout << "Input VPD commit time: " << (double)(clock() - commit_t) / CLOCKS_PER_SEC << endl;
	
	return make_pair(r_f, r_f2);
}

bool check_commit(Ec1 digest, Ec1 digesta, Ec1 digest2, Ec1 digest2a){
	Fp12 ea1, ea2;
	
	opt_atePairing(ea1, g2, digesta);
	opt_atePairing(ea2, g2a, digest);

	Fp12 ea3, ea4;
	
	opt_atePairing(ea3, g2, digest2a);
	opt_atePairing(ea4, g2a, digest2);
	
	return (ea1 == ea2) && (ea3 == ea4);
}

void prove(vector<mpz_class> r, mpz_class& ans, vector<mpz_class>& input, vector<mpz_class> &input2, vector<Ec1>& witness, vector<Ec1>& witnessa, mpz_class r_f, mpz_class r_f2, mpz_class Z){
	vector<mpz_class> coeffs = input;
	coeffs[0] = (coeffs[0] + (Z * input2[0]) % p) % p;
	coeffs[1 << (NumOfVar - 1)] = (coeffs[1 << (NumOfVar - 1)] + (Z * input2[1]) % p) % p;

	clock_t prove_t = clock();

	std::vector<mpz_class> t(NumOfVar);
	for(int i = 0; i < t.size(); i++)
		mpz_urandomm(t[i].get_mpz_t(), r_state, p.get_mpz_t());
	
	assert(coeffs.size() >= 1);
	coeffs[coeffs.size() - 1] = ((r_f + Z * r_f2) % p);


	for(int i = 0; i < coeffs.size(); i++)
		if(coeffs[i] < 0)
			coeffs[i] += p;
			
	std::vector<mpz_class> ans_pre;
	ans_pre.resize((int)pow(2, NumOfVar));
	ans_pre[0] = 1;
	ans = coeffs[0];
	for(int i = 0; i < NumOfVar; i++){
		for(int j = (int)pow(2, i); j < (int)pow(2, i+1); j++){
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
	for(int i = 0; i < NumOfVar; i++){
		for(int j = 0; j < pow(2, r.size() - i - 1); j++){
			witness_coeffs[start_index + j] = temp_coeffs[pow(2, r.size() - i - 1) + j] % p;
			temp_coeffs[j] = (temp_coeffs[j] + (temp_coeffs[pow(2, r.size() - i - 1) + j] * r[NumOfVar - 1 - i]) % p) % p;

		}
		start_index += pow(2, r.size() - i - 1);
	}
	for(int i = 0; i < witness_coeffs.size(); i++)
		if(witness_coeffs[i] < 0)
			witness_coeffs[i] += p;

	const mie::Vuint tempa(a.get_str().c_str());

	std::vector<mpz_class> scalar_pow;
	scalar_pow.resize(multi_scalar_w);
	for(int i = 0; i < (1 << (NumOfVar - 1)) / multi_scalar_w + 1; ++i)
	{
		int temp = 0;
		for(int k = NumOfVar - 1; k >= 0; --k)
		{
			if(i >= ((1 << k) / multi_scalar_w + 1))
			{
				break;
			}
			for(int j = 0; j < multi_scalar_w; ++j)
			{
				int id = i * multi_scalar_w + j;
				if(id >= (1 << k))
				{
					scalar_pow[j] = 0;
				}
				else
				{
					scalar_pow[j] = witness_coeffs[id + temp];
				}
			}
			temp += (1 << k);
			witness[k] += multi_scalar_calc(i, (1 << k), scalar_pow);
		}
	}

	for(int k = NumOfVar - 1; k >= 0; --k)
	{
		mie::Vuint temptk(t[k].get_str().c_str());
		witness[k] += pub_g1[1 << NumOfVar] * temptk;
		witnessa[k] = witness[k] * tempa;
	}


	mpz_class tmp = coeffs[1 << NumOfVar];
	for(int i = 0; i < NumOfVar; i++)
		tmp += (t[i] * r[i]) % p;
	witness[NumOfVar] = pre_exp(g1_pre, tmp);
	clock_t zkt = clock();
	for(int i = 0; i < NumOfVar; i++){
		mie::Vuint tempti(t[i].get_str().c_str());
		witness[NumOfVar] -= pub_g1[1 << i] * tempti;
	}
	witnessa[NumOfVar] = witness[NumOfVar] * tempa;

	cout << "Input VPD prove time: " << (double)(clock() - prove_t) / CLOCKS_PER_SEC << endl;	
}

bool verify(vector<mpz_class> r, Ec1 digest, Ec1 digest2, mpz_class Z, mpz_class& ans, vector<Ec1>& witness, vector<Ec1>& witnessa){
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

	mie::Vuint tempz(Z.get_str().c_str());
	Ec1 temp2 = pre_exp(g1_pre, ans);

	opt_atePairing(ea3, g2, digest + digest2 * tempz - temp2);

	for(int i = 0; i < r.size() + 1; i++){
		if(i == r.size()){
			opt_atePairing(temp[i], pub_g2[i], witness[i]);
		}
		if(i < r.size()){
			Ec2 temp3 = pub_g2[i] - pre_exp(g2_pre, r[i]);
		
			opt_atePairing(temp[i], temp3, witness[i]);
		}
		ea4 *= temp[i];
	}





	if(ea3 != ea4) {
		cout << "final error" << endl;
		flag = 0;
	}

	cout << "Input VPD verify time: "<<(double)(clock() - verify_t) / CLOCKS_PER_SEC << endl;
	return flag;
	
}

void environment_init()
{
	p.set_str("16798108731015832284940804142231733909759579603404752749028378864165570215949",10);
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
}

}
