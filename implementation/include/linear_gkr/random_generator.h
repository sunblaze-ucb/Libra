#ifndef __random_gen
#define __random_gen

#include <cstdlib>
#include <ctime>
#include <cstring>
#include <cassert>
/*
This is a implementation of random oracle
*/

class random_oracle
{
private:
public:
	random_oracle(){}
	~random_oracle(){}
	unsigned long long query(unsigned long long x)
	{
		//most simple impl, debug purpose only
		//it should be a hash function
		#ifdef __debug
		srand(x);
		return rand();
		#endif
		//if __debug is not define and no secure implemetation, the program will throw a exception
		assert(false);
	}
};
#endif