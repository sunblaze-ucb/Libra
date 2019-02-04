#include <cstdio>
#include "traditional_gkr/verifier_traditional.h"
#include "traditional_gkr/prover_clogc.h"
#include "linear_gkr/prime_field.h"

verifier* v;
prover* p;


int input_size;
int n_copy, log_n_copy;

void run_verifications(const char *input_path)
{
	bool final_result = true;
	double total_time = 0;
	p = new prover[n_copy];
	v = new verifier[n_copy];
	for(int i = 0; i < n_copy; ++i)
	{
		p[i].total_time = 0;
		v[i].get_prover(&p[i]);
		FILE *input_file = fopen(input_path, "r");
		v[i].read_circuit_from_FILE(input_file);
		fclose(input_file);
		p[i].get_circuit(v[i].C);
		bool result = v[i].verify();
		final_result &= result;
		total_time += p[i].total_time;
	}
	if(final_result)
	{
		printf("Verification Pass, total time %f\n", (float)total_time);
	}
	else
	{
		printf("Verification Fail.\n");
	}
	delete[] p;
	delete[] v;
}

int main(int argc, char** argv)
{
	sscanf(argv[2], "%d", &n_copy);
	sscanf(argv[3], "%d", &log_n_copy);
	prime_field::init("16798108731015832284940804142231733909759579603404752749028378864165570215949", 10);
	run_verifications(argv[1]);
	return 0;
}
