#include <cstdio>

int main()
{
	int d = 3;
	int w[] = {20, 20, 19};
	printf("%d\n", d);
	int n = (1 << w[0]);
	printf("%d ", n);
	for(int j = 0; j < n; ++j)
	{
		printf("%d %d %d %d ", 3, j, j + 1, 0);
	}
	printf("\n");
	n = (1 << w[1]);
	printf("%d ", n);
	for(int j = 0; j < n; ++j)
	{
		printf("%d %d %d %d ", 0, j, j, j);
	}
	printf("\n");
	n = (1 << w[2]);
	printf("%d ", n);
	for(int j = 0; j < n; ++j)
	{
		printf("%d %d %d %d ", 0, j, (j << 1), (j << 1 | 1));
	}
	printf("\n");

	return 0;
}
