#include <cstdio>
#include <cstdlib>
#include <utility>
#include <cassert>
#include <vector>

using namespace std;

int main(int argc, char **argv)
{
	int mode;
	sscanf(argv[1], "%d", &mode);
	if(mode == 0) //three layers of simple circuit
	{
		int d, w;
		sscanf(argv[2], "%d", &d);
		sscanf(argv[3], "%d", &w);
		printf("%d\n", d);
		int n = (1 << w);
		printf("%d ", n);
		for(int j = 0; j < n; ++j)
		{
			printf("%d %d %d %d ", 3, j, rand(), 0);
		}
		printf("\n");
		for(int k = 1; k < d - 1; ++k)
		{
			n = (1 << w);
			printf("%d ", n);
			for(int j = 0; j < n; ++j)
			{
				printf("%d %d %d %d ", 8, j, j, 0);
			}
			printf("\n");
		}
		n = (1 << (w - 1));
		printf("%d ", n);
		for(int j = 0; j < n; ++j)
		{
			printf("%d %d %d %d ", 10, j, (j << 1), 0);
		}
		printf("\n");
	}
	else if(mode == 1) //blocked circuit
	{
		FILE *meta;
		meta = fopen("meta_data.txt", "w");
		int d, w;
		sscanf(argv[2], "%d", &d);
		sscanf(argv[3], "%d", &w);
		assert(w > 5);
		printf("%d\n", d);
		int n = (1 << w);
		int block_number = (1 << (w - 5));
		assert(n >= block_number);
		int block_size = n / block_number;
		int str_length;
		str_length = 0;
		str_length += printf("%d ", n);
		for(int j = 0; j < n; ++j)
		{
			str_length += printf("%d %d %010d %d ", 3, j, rand(), 0);
		}
		str_length += printf("\n");
		fprintf(meta, "%d %d %d %d %d\n", block_size, block_number, 5, w - 5, str_length);
		str_length = 0;
		for(int k = 1; k < d - 1; ++k)
		{
			n = (1 << w);
			str_length += printf("%d ", n);
			for(int j = 0; j < n; ++j)
			{
				str_length += printf("%d %d %d %d ", 0, j, j, j);
			}
			str_length += printf("\n");
			block_size = n / block_number;
			assert(n >= block_size);
			fprintf(meta, "%d %d %d %d %d\n", block_size, block_number, 5, w - 5, str_length);
			str_length = 0;
		}
		n = (1 << (w - 1));
		block_size = n / block_number;
		str_length += printf("%d ", n);
		for(int i = 0; i < block_number; ++i)
		{
			for(int j = 0; j < block_size; ++j)
			{
				int id = (i * block_size) + j;
				int nxta = (i * block_size * 2) + j * 2;
				int nxtb = (i * block_size * 2) + j * 2 + 1;
				str_length += printf("%d %d %d %d ", 1, id, nxta, nxtb);
			}
		}
		str_length += printf("\n");
		assert(n >= block_size);
		fprintf(meta, "%d %d %d %d %d\n", block_size, block_number, 4, w - 5, str_length);
		str_length = 0;
		fclose(meta);
	}
	else if(mode == 2) //matrix mul with addition tree
	{
		FILE *meta;
		meta = fopen("meta_data.txt", "w");

		int mat_sz;
		sscanf(argv[2], "%d", &mat_sz);

		assert(__builtin_popcount(mat_sz) == 1);

		int log_mat_sz = 0;
		while(mat_sz != (1 << log_mat_sz))
			log_mat_sz++;
		
		//input layer
		int block_number = mat_sz * mat_sz;
		int block_size;
		vector<vector<int> > A, B;
		A.resize(mat_sz), B.resize(mat_sz);

		for(int i = 0; i < mat_sz; ++i)
		{
			A[i].resize(mat_sz);
			B[i].resize(mat_sz);
			for(int j = 0; j < mat_sz; ++j)
				A[i][j] = rand() % 10, B[i][j] = rand() % 10;
		}
		//expand the input for trivial data parallel

		//input layer
		int str_length = 0;
		printf("%d\n", 1 + 1 + log_mat_sz);
		str_length += printf("%d ", mat_sz * mat_sz * mat_sz * 2);
		for(int i = 0; i < mat_sz; ++i)
		{
			for(int j = 0; j < mat_sz; ++j)
			{
				for(int k = 0; k < mat_sz; ++k)
				{
					str_length += printf("%d %d %010d %d ", 3, 2 * (i * mat_sz * mat_sz + j * mat_sz + k) + 0, A[i][k], 0);
					str_length += printf("%d %d %010d %d ", 3, 2 * (i * mat_sz * mat_sz + j * mat_sz + k) + 1, B[k][j], 0);
				}
			}
		}
		str_length += printf("\n");
		fprintf(meta, "%d %d %d %d %d\n", mat_sz * 2, block_number, log_mat_sz + 1, 2 * log_mat_sz, str_length);

		//mult
		str_length = printf("%d ", mat_sz * mat_sz * mat_sz);
		for(int i = 0; i < mat_sz; ++i)
		{
			for(int j = 0; j < mat_sz; ++j)
			{
				for(int k = 0; k < mat_sz; ++k)
				{
					int id = i * mat_sz * mat_sz + j * mat_sz + k;
					str_length += printf("%d %d %d %d ", 1, id, id * 2, id * 2 + 1);
				}
			}
		}
		str_length += printf("\n");

		fprintf(meta, "%d %d %d %d %d\n", mat_sz, block_number, log_mat_sz, 2 * log_mat_sz, str_length);

		//addition tree
		int num_leaves = mat_sz;
		for(int dep = 0; dep < log_mat_sz; ++dep)
		{
			num_leaves /= 2;
			str_length = printf("%d ", num_leaves * mat_sz * mat_sz);
			for(int i = 0; i < mat_sz; ++i)
			{
				for(int j = 0; j < mat_sz; ++j)
				{
					for(int k = 0; k < num_leaves; ++k)
					{
						int id = i * mat_sz * num_leaves + j * num_leaves + k;
						str_length += printf("%d %d %d %d ", 0, id, id * 2, id * 2 + 1);
					}
				}
			}
			str_length += printf("\n");
			fprintf(meta, "%d %d %d %d %d\n", num_leaves, block_number, log_mat_sz - dep - 1, 2 * log_mat_sz, str_length);
		}
		fclose(meta);
	}
	else if(mode == 3) //matrix mul with non-expanded input
	{
		
	}
	else if(mode == 4) //matrix mul with both non-expanded input and summation gate
	{
		int mat_sz;
		sscanf(argv[2], "%d", &mat_sz);

		assert(__builtin_popcount(mat_sz) == 1);

		int log_mat_sz = 0;
		while(mat_sz != (1 << log_mat_sz))
			log_mat_sz++;
		
		//input layer
		int block_number = mat_sz * mat_sz;
		int block_size;
		vector<vector<int> > A, B;
		A.resize(mat_sz), B.resize(mat_sz);

		for(int i = 0; i < mat_sz; ++i)
		{
			A[i].resize(mat_sz);
			B[i].resize(mat_sz);
			for(int j = 0; j < mat_sz; ++j)
				A[i][j] = rand() % 10, B[i][j] = rand() % 10;
		}
		//input layer
		printf("%d\n", 1 + 1 + 1);
		printf("%d ", mat_sz * mat_sz * 2);
		for(int i = 0; i < mat_sz; ++i)
		{
			for(int j = 0; j < mat_sz; ++j)
			{
				printf("%d %d %010d %d ", 3, i * mat_sz + j, A[i][j], 0);
			}
		}
		for(int i = 0; i < mat_sz; ++i)
		{
			for(int j = 0; j < mat_sz; ++j)
			{
				printf("%d %d %010d %d ", 3, mat_sz * mat_sz + i * mat_sz + j, B[i][j], 0);
			}
		}
		printf("\n");

		//mult
		printf("%d ", mat_sz * mat_sz * mat_sz);
		for(int i = 0; i < mat_sz; ++i)
		{
			for(int j = 0; j < mat_sz; ++j)
			{
				for(int k = 0; k < mat_sz; ++k)
				{
					int id = i * mat_sz * mat_sz + j * mat_sz + k;
					int a = i * k, b = k * j;
					printf("%d %d %d %d ", 1, id, a, b);
				}
			}
		}
		printf("\n");

		

		//addition tree
		printf("%d ", mat_sz * mat_sz);
		for(int i = 0; i < mat_sz; ++i)
		{
			for(int j = 0; j < mat_sz; ++j)
			{
				printf("5 %d %d %d ", i * mat_sz + j, i * mat_sz * mat_sz + j * mat_sz, i * mat_sz * mat_sz + j * mat_sz + mat_sz);
			}
		}
	}
	else if(mode == 5) //addition tree test
	{
		int d, w;
		sscanf(argv[2], "%d", &d);
		sscanf(argv[3], "%d", &w);
		printf("%d\n", d);
		int n = (1 << w);
		printf("%d ", n);
		for(int j = 0; j < n; ++j)
		{
			printf("%d %d %d %d ", 3, j, 1, 0);
		}
		printf("\n");
		printf("2 5 0 0 0 5 1 0 0\n");
	}
	return 0;
}
