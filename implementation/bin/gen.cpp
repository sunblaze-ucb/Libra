#include <cstdio>
#include <cassert>

int main(int argc, char **argv)
{
	int mode;
	sscanf(argv[1], "%d", &mode);
	if(mode == 0)
	{
		int d, w;
		sscanf(argv[2], "%d", &d);
		sscanf(argv[3], "%d", &w);
		printf("%d\n", d);
		int n = (1 << w);
		printf("%d ", n);
		for(int j = 0; j < n; ++j)
		{
			printf("%d %d %d %d ", 3, j, j + 1, 0);
		}
		printf("\n");
		for(int k = 1; k < d - 1; ++k)
		{
			n = (1 << w);
			printf("%d ", n);
			for(int j = 0; j < n; ++j)
			{
				printf("%d %d %d %d ", 0, j, j, j);
			}
			printf("\n");
		}
		n = (1 << (w - 1));
		printf("%d ", n);
		for(int j = 0; j < n; ++j)
		{
			printf("%d %d %d %d ", 1, j, (j << 1), (j << 1 | 1));
		}
		printf("\n");
	}
	else
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
			str_length += printf("%d %d %010d %d ", 3, j, j + 1, 0);
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
	return 0;
}
