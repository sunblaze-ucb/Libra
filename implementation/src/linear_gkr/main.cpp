#define __debug__
#include "linear_gkr/verifier.h"
#include "linear_gkr/prover.h"
#include "linear_gkr/prime_field.h"

verifier v;
prover p;

int main()
{
	prime_field::field_element test;
	test.set_value(1);
	test = test + test;
	prime_field::init("16798108731015832284940804142231733909759579603404752749028378864165570215949", 10);
	v.get_prover(&p);
	return 0;
}