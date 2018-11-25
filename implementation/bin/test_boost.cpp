#include <boost/multiprecision/cpp_int.hpp>
#include <boost/random.hpp>

int main()
{
	using namespace boost::multiprecision;
	using namespace boost::random;

	independent_bits_engine<mt19937, 256, cpp_int> gen;

	std::cout << gen() << std::endl;
	return 0;
}