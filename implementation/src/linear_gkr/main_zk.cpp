#define __debug__
#define __timer__
#include "linear_gkr/zk_verifier.h"
#include "linear_gkr/zk_prover.h"
#include "linear_gkr/prime_field.h"
#include <iostream>
zk_verifier v;
zk_prover p;

int main()
{
	//std::cout << "hello world" << std::endl;

	prime_field::init("16798108731015832284940804142231733909759579603404752749028378864165570215949", 10);
	p.total_time = 0;
	v.get_prover(&p);
	//std::cout << "come in" << std::endl;
	v.read_circuit("test_circuit.txt");
	//std::cout << "after readfile" << std::endl;
	p.get_circuit(v.C);
	bool result = v.verify();
	printf("%s\n", result ? "Pass" : "Fail");
	return 0;
}