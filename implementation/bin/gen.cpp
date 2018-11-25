#include <cstdio>

int main(int argc, char **argv)
{
	int d, w;
	sscanf(argv[1], "%d", &d);
	sscanf(argv[2], "%d", &w);
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
		printf("%d %d %d %d ", 0, j, (j << 1), (j << 1 | 1));
	}
	printf("\n");

	return 0;
}
