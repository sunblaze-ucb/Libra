#include <cstdio>
#include "traditional_gkr/verifier_traditional.h"
#include "traditional_gkr/prover_clogc.h"
#include "linear_gkr/prime_field.h"

verifier* v;
prover* p;

int num_blocks;

char** inputs;
char** input_ptr;
int input_size;
void run_verifications()
{
	bool final_result = true;
	double total_time = 0;
	p = new prover[num_blocks];
	v = new verifier[num_blocks];
	for(int i = 0; i < num_blocks; ++i)
	{
		p[i].total_time = 0;
		v[i].get_prover(&p[i]);
		v[i].read_circuit_from_string(inputs[i]);
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

void redistribute_circuit(const char *path)
{
	FILE *circuit, *meta;
	circuit = fopen(path, "r");
	meta = fopen("meta_data.txt", "r");

	int d;
	fscanf(circuit, "%d", &d);
	input_size = 0;
	int block_size, block_number, binary_bs, binary_bn;
	int str_length;
	for(int i = 0; i < d; ++i)
	{
		fscanf(meta, "%d%d%d%d%d", &block_size, &block_number, &binary_bs, &binary_bn, &str_length);
		input_size += str_length / block_number * 2;
	}
	rewind(meta);
	inputs = new char*[block_number];
	input_ptr = new char*[block_number];
	for(int i = 0; i < block_number; ++i)
		inputs[i] = new char[input_size * 2];
	for(int i = 0; i < block_number; ++i)
		input_ptr[i] = inputs[i];
	for(int i = 0; i < block_number; ++i)
	{
		input_ptr[i] += sprintf(input_ptr[i], "%d\n", d);
	}

	for(int i = 0; i < d; ++i)
	{
		int n;
		fscanf(circuit, "%d", &n);
		int p_bs = block_size, p_bn = block_number;
		fscanf(meta, "%d%d%d%d%d", &block_size, &block_number, &binary_bs, &binary_bn, &str_length);
		num_blocks = block_number;
		assert(block_number == (1 << binary_bn));
		assert(block_size == (1 << binary_bs));
		assert(n == block_number * block_size);
		for(int j = 0; j < block_number; ++j)
		{
			input_ptr[j] += sprintf(input_ptr[j], "%d ", block_size);
			int offset = j * block_size;
			int p_offset = j * p_bs;
			for(int k = 0; k < block_size; ++k)
			{
				int ty, id, nxta, nxtb;
				fscanf(circuit, "%d%d%d%d", &ty, &id, &nxta, &nxtb);
				switch(ty)
				{
					case 0:
						id -= offset;
						nxta -= p_offset;
						nxtb -= p_offset;
						break;
					case 1:
						id -= offset;
						nxta -= p_offset;
						nxtb -= p_offset;
						break;
					case 2:
						break;
					case 3:
						id -= offset;
						break;
				}
				input_ptr[j] += sprintf(input_ptr[j], "%d %d %d %d ", ty, id, nxta, nxtb);
			}
			input_ptr[j] += sprintf(input_ptr[j], "\n");
		}
	}
	fclose(circuit);
	fclose(meta);
}

int main(int argc, char** argv)
{
	redistribute_circuit("test_circuit.txt");
	prime_field::init("16798108731015832284940804142231733909759579603404752749028378864165570215949", 10);
	run_verifications();
	for(int i = 0; i < num_blocks; ++i)
	{
		delete[] inputs[i];
	}
	delete[] inputs;
	delete[] input_ptr;
	return 0;
}
