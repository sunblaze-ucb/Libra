#include <gmp.h>
#include <gmpxx.h>

int main()
{
	mpz_class* a = new mpz_class[10];
	for(int i = 0; i < 10; ++i)
		a[i] = mpz_class(i);
	delete[] a;
	return 0;
}