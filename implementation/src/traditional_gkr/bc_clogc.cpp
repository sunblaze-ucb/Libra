#include <cstdio>
#include "traditional_gkr/verifier_bc_clogc.h"
#include "traditional_gkr/prover_bc_clogc.h"
#include "linear_gkr/prime_field.h"

verifier v;
prover p;

int input_size;
int n_copy, log_n_copy;
void run_verifications(const char *path)
{
	bool final_result = true;
	double total_time = 0;
	v.set_blocks(n_copy, log_n_copy);
	v.get_prover(&p);
	p.total_time = 0;
	v.read_circuit(path, n_copy);
	p.get_circuit(v.C);
	final_result &= v.verify();
	total_time = p.total_time;
	if(final_result)
	{
		printf("Verification Pass, total time %f\n", (float)total_time);
	}
	else
	{
		printf("Verification Fail.\n");
	}
}

int main(int argc, char** argv)
{
	sscanf(argv[2], "%d", &n_copy);
	sscanf(argv[3], "%d", &log_n_copy);
	prime_field::init("16798108731015832284940804142231733909759579603404752749028378864165570215949", 10);
	run_verifications(argv[1]);
	return 0;
}
