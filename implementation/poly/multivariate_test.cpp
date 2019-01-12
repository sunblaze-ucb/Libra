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

#define P 254

Ec1 g1;
Ec2 g2, ga;
mpz_class p;
unsigned long int seed;
gmp_randstate_t r_state;

int W = 2;
int D;

vector<vector<Ec1> > pubs_g1,pubs_g1a;

vector<vector<vector<Ec1> > > pubs_g1_pre, pubs_g1a_pre;
vector<vector<Ec1> > pubs_g1_pre2, pubs_g1a_pre2;

vector<Ec2> pubs_g2;

vector<Ec1> g1_pre, g1a_pre;
vector<Ec2> g2_pre;

void precompute_g1(mpz_class a){
	g1_pre.resize(P);
	g1a_pre.resize(P);
	g2_pre.resize(P);
	g1_pre[0] = g1;
	g2_pre[0] = g2;
	const mie::Vuint temp(a.get_str().c_str());
	g1a_pre[0] = g1*temp;
	for(int i=1;i<P;i++){
		g1_pre[i]=g1_pre[i-1]+g1_pre[i-1];
		g2_pre[i] = g2_pre[i-1]+g2_pre[i-1];
		//g1a_pre[i]=g1a_pre[i-1]+g1a_pre[i-1];
	}

}

void precompute_pub(){
	pubs_g1_pre.resize(D+1);
	pubs_g1a_pre.resize(D+1);
	for(int i=0;i<1;i++){
		pubs_g1_pre[i].resize((int)pow(2,i));
		pubs_g1a_pre[i].resize((int)pow(2,i));
		for(int j=0;j<(int)pow(2,i);j++){
			pubs_g1_pre[i][j].resize(P);
			pubs_g1a_pre[i][j].resize(P);	
		}
	
	}
	
	for(int i=0;i<1;i++){
		for(int j=0;j<(int)pow(2,i);j++){
			pubs_g1_pre[i][j][0] = pubs_g1[i][j];
			pubs_g1a_pre[i][j][0] = pubs_g1a[i][j];	
			for(int k=1;k<P;k++){
				pubs_g1_pre[i][j][k] = pubs_g1_pre[i][j][k-1]*2;
				pubs_g1a_pre[i][j][k] = pubs_g1a_pre[i][j][k-1]*2;
			}
		}
	
	}
	

	
	return;

}


void precompute_pub2(){
	pubs_g1_pre2.resize(D+1);
	pubs_g1a_pre2.resize(D+1);
	for(int i=1;i<D+1;i++){
		pubs_g1_pre2[i].resize((int)pow(2,i)/W);
		pubs_g1a_pre2[i].resize((int)pow(2,i)/W);
	
	}
	

	
	for(int i=1;i<D+1;i++){
		for(int j=0;j<(int)pow(2,i)/W;j++){
			vector<Ec1> temp(W);
			for(int k=0;k<W;k++)
				temp[k] = pubs_g1[i][j*W+k];	
			pubs_g1_pre2[i][j]=temp[0]+temp[1];

			
		}
	
	}
	
	for(int i=1;i<D+1;i++){
		for(int j=0;j<(int)pow(2,i)/W;j++){
			vector<Ec1> temp(W);
			for(int k=0;k<W;k++)
				temp[k] = pubs_g1a[i][j*W+k];
			pubs_g1a_pre2[i][j]=temp[0]+temp[1];
				
		}
	
	}
	
	
	
	return;

}




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



template< class T >
T pre_exp(vector<T>& pre, mpz_class n){
	T temp = pre[0]*0;
	int length = mpz_sizeinbase(n.get_mpz_t(), 2);
	//cout << "length = " << length << endl;
	for(int i=0;i<length;i++){
		if(mpz_tstbit(n.get_mpz_t(),i)==1)
			temp = temp + pre[i];		
	}
	
	
	return temp;

}


void keygen(int d){
	clock_t t1 = clock();
	D=d;
	
	//secret keys
	vector<mpz_class> s(D);
	
	for(int i=0;i<D;i++)
		mpz_urandomm(s[i].get_mpz_t(),r_state,p.get_mpz_t());
	
	mpz_class a;
	
	mpz_urandomm(a.get_mpz_t(),r_state,p.get_mpz_t());
	
	
	//compute public keys
	precompute_g1(a);
	
	
	pubs_g1.resize(D+1);
	pubs_g1a.resize(D+1);
	pubs_g2.resize(D);

	
	for(int i=0;i<D+1;i++){
		pubs_g1[i].resize((int)pow(2,i));
		pubs_g1a[i].resize((int)pow(2,i));
	}
	
	pubs_g1[0][0] = g1;
	pubs_g1a[0][0] = g1a_pre[0];
	
	
	
	{
		const mie::Vuint temp((a.get_str()).c_str());
		ga = g2*temp;
	}
	
	
	vector<vector<mpz_class> > vars(D+1);
	for(int i=0;i<D+1;i++){
		vars[i].resize((int)pow(2,i));
		
	}
	
	vars[0][0]=1;
	
	for(int i=1;i<D+1;i++){
		for(int j=0;j<(int)pow(2,i-1);j++){
			vars[i][2*j+1] = (vars[i-1][j]*s[D-i])%p;
			
			vars[i][2*j] = (vars[i-1][j]*(1-s[D-i]))%p;
			
			//std::cout << "vars[i][2*j+1] = " << (vars[i][2*j+1].get_str()).c_str() << std::endl;
		}
	}
	
	
	for(int i=1;i<D+1;i++){
		for(int j=0;j<(int)pow(2,i-1);j++){
			
			if(vars[i][2*j+1]<0)
				vars[i][2*j+1]+=p;

				
			pubs_g1[i][2*j+1] = pre_exp(g1_pre,vars[i][2*j+1]);
			pubs_g1[i][2*j] = pubs_g1[i-1][j]-pubs_g1[i][2*j+1];
			pubs_g1a[i][2*j+1] = pre_exp(g1_pre,(vars[i][2*j+1]*a)%p);
			pubs_g1a[i][2*j] = pubs_g1a[i-1][j]-pubs_g1a[i][2*j+1];
		}
	}
	
	

	for(int i=0;i<D;i++){
		const mie::Vuint temp((s[i].get_str()).c_str());
		pubs_g2[i]=g2*temp;
	}
	
	vars.resize(0);
	
	precompute_pub();
	precompute_pub2();
	
	cout<<"pub time: "<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s\n";
	
	return;
}

void commit(Ec1& digest, Ec1& digesta, vector<mpz_class>& input){
	vector<mpz_class> coeffs=input;
	
	int d = ceil(log2(input.size()));
	
	clock_t t1 = clock();
	
	//compute digest pub
	
	for(int i=0;i<coeffs.size();i++)
		if(coeffs[i]<0)
			coeffs[i]+=p;
	digest = multi_scalar(pubs_g1[d], pubs_g1_pre2[d], coeffs);
	digesta = multi_scalar(pubs_g1a[d], pubs_g1a_pre2[d], coeffs);
	
	
	cout<<"digest time: "<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s\n";
	
	return;
	
}


bool check_commit(Ec1 digest, Ec1 digesta){
	Fp12 ea1, ea2;
	
	opt_atePairing(ea1, g2, digesta);
	opt_atePairing(ea2, ga, digest);
	
	
	return (ea1==ea2);
	
}

void prove(vector<mpz_class> r, mpz_class& ans, vector<mpz_class>& input, vector<Ec1>& witness, vector<Ec1>& witnessa){
	
	
	//answer
	vector<mpz_class> products((int)pow(2,r.size())), coeffs = input;
	products[0] = 1;
	int index = 1;
	for(int i=0;i<r.size();i++){
		for(int j=0;j<pow(2,i);j++){
			products[index+j] = (products[j]*r[i])%p;
			products[j] = (products[j]*(1-r[i]))%p;
		}
		index+=(int)pow(2,i);
	}

	ans = 0;
	for(int i=0;i<coeffs.size();i++){
		ans+=products[i]*coeffs[i];
	}


	ans%=p;
	if(ans<0)
		ans+=p;
	
	//witness coeffs
	
	clock_t t1 = clock();
	vector<mpz_class> witness_coeffs((int)pow(2,r.size())), temp_coeffs = coeffs;
	coeffs.resize(0);
	temp_coeffs.resize((int)pow(2,r.size()));
	
	int start_index = 0;
	
	for(int i=0;i<r.size();i++){
		for(int j=0;j<pow(2,r.size()-i-1);j++){
			witness_coeffs[start_index+j] = (-temp_coeffs[2*j]+temp_coeffs[2*j+1])%p;
			temp_coeffs[j] = (-temp_coeffs[2*j]*(r[i]-1)+temp_coeffs[2*j+1]*r[i])%p;
		}
		temp_coeffs.resize(pow(2,r.size()-i-1));
		start_index+=pow(2,r.size()-i-1);
	}
	

	
	for(int i=0;i<witness_coeffs.size();i++){
		if(witness_coeffs[i]<0)
			witness_coeffs[i]+=p;
		
	}
	
	
	//witness
	witness.resize(r.size()), witnessa.resize(r.size());
	
	start_index = 0;
	
	
	for(int i=0;i<r.size()-1;i++){
		vector<mpz_class> temp_e((int)pow(2,r.size()-i-1));
		for(int j=0;j<temp_e.size();j++){
			temp_e[j] = witness_coeffs[start_index+j];
		}
		
		witness[i] = multi_scalar(pubs_g1[r.size()-i-1], pubs_g1_pre2[r.size()-i-1], temp_e);
		if(i!=0)
			witnessa[i] = multi_scalar(pubs_g1a[r.size()-i-1], pubs_g1a_pre2[r.size()-i-1], temp_e);
		
		start_index+=pow(2,r.size()-i-1);
		
	}
	
	
	
	for(int i=r.size()-1;i<r.size();i++){
		witness[i] = g1*0;
		witnessa[i] = g1*0;
		for(int j=0;j<pow(2,r.size()-i-1);j++){
			if(witness_coeffs[start_index+j]!=0){
			
				mpz_class temp1 = witness_coeffs[start_index+j]; 
				mpz_class temp2 = temp1-p;
				if(temp1.get_str().length()<temp2.get_str().length()){
					//const mie::Vuint temp((temp1.get_str()).c_str());
					witness[i] += pre_exp(pubs_g1_pre[r.size()-i-1][j],temp1);	
					witnessa[i] += pre_exp(pubs_g1a_pre[r.size()-i-1][j],temp1);
					
				}
				else{
					temp2 = (-temp2);
					witness[i] -= pre_exp(pubs_g1_pre[r.size()-i-1][j],temp2);	
					witnessa[i] -= pre_exp(pubs_g1a_pre[r.size()-i-1][j],temp2);
						

				
				}
			}
		}
		
		start_index+=pow(2,r.size()-i-1);
		
	}
	
	

	cout<<"prover time: "<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s\n";
	
	return;
}

bool verify(vector<mpz_class> r, Ec1 digest, mpz_class& ans, vector<Ec1>& witness, vector<Ec1>& witnessa){
	clock_t t1 = clock();
	
	
	Fp12 ea1, ea2;
	
	bool flag=1;
	
	for(int i=0;i<r.size();i++){
		if(i!=0){
			opt_atePairing(ea1, g2, witnessa[i]);
			opt_atePairing(ea2, ga, witness[i]);
			
			
			//cout<<(ea1==ea2)<<endl;
			
			flag = (flag&(ea1==ea2));
		}
		
		
	}
	
	//cout<<"flag: "<<flag<<endl;
	
	Fp12 e1,e3=1;
	vector<Fp12> e2(r.size());
	//const mie::Vuint templ(ans.get_str().c_str());
	
	Ec1 temp2 = pre_exp(g1_pre, ans);
	
	
	opt_atePairing(e1,g2, digest-temp2);

	for(int i=0;i<r.size();i++){
		//const mie::Vuint temp(r[i].get_str().c_str());
	
		
		Ec2 temp1 = pubs_g2[D-r.size()+i]-pre_exp(g2_pre,r[i]);
		
		opt_atePairing(e2[i],temp1, witness[i]);
		
		e3*=e2[i];
		
	}

	//cout<<"verification:"<<(e1==e3)<<"\n";	
	
	cout<<"verify time: "<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s\n";
	//vtime+=(double)(clock()-t1)/CLOCKS_PER_SEC;
	return (flag&(e1==e3));
	
}


int main(int argc, char** argv){
	
	
	
	seed = rand();
    gmp_randinit_default (r_state);
    gmp_randseed_ui(r_state, seed);
	p.set_str("16798108731015832284940804142231733909759579603404752749028378864165570215949",10);
	
	
	//bilinear g1 g2
	bn::CurveParam cp = bn::CurveFp254BNb;
	Param::init(cp);
	const Point& pt = selectPoint(cp);
	g2=Ec2(
		Fp2(Fp(pt.g2.aa), Fp(pt.g2.ab)),
		Fp2(Fp(pt.g2.ba), Fp(pt.g2.bb))
	);
	g1=Ec1(pt.g1.a, pt.g1.b);
	
	
	
	int d = atoi(argv[1]);
	
	keygen(d);
	
	int N = (int)pow(2,D);
	vector<mpz_class> input(N);
	for(int i=0;i<input.size();i++){
		input[i] = rand();
	}
	
	
	
	Ec1 digest, digesta;
	
	commit(digest,digesta,input);
	cout<<"check commit: "<<check_commit(digest,digesta)<<endl;
	
	vector<Ec1> proof, proofa;
	vector<mpz_class> r(d);
	mpz_class ans;
	
	for(int i=0;i<r.size();i++)
		mpz_urandomm(r[i].get_mpz_t(),r_state,p.get_mpz_t());
	
	prove(r, ans, input, proof, proofa);
	
	bool tf = verify(r,digest,ans,proof,proofa);
	cout<<"verify: "<<tf<<endl;
	
	
	return 1;
	
}