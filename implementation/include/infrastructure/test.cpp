
#include "RS_polynomial.h"
#include <utility>
#include <chrono>
#include <iostream>
const int lg_order = 10;
const int lg_coef = 5;
const int num_coef = 1 << lg_coef;
const int num_order = 1 << lg_order;
prime_field::field_element coef[num_coef << 1];

int main()
{
    for(int i = 0; i < (num_coef << 1); ++i)
    {
        coef[i] = prime_field::random();
    }
    auto rou = prime_field::get_root_of_unity(lg_order);
    
	std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    auto res = fast_fourier_transform(coef, num_coef, num_order, rou);
    auto coe_res = inverse_fast_fourier_transform(res, num_coef, num_order, rou);
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	std::cerr << "total evaluation time: " << time_span.count() << " seconds." << std::endl;
    
    auto x = prime_field::field_element(1);
    for(int i = 0; i < num_order; ++i)
    {
        auto sum = prime_field::field_element(0);
        auto cur_x = prime_field::field_element(1);
        for(int j = 0; j < num_coef; ++j)
        {
            sum = sum + coef[j] * cur_x;
            cur_x = cur_x * x;
        }
        x = x * rou;
        assert(sum == res[i]);
    }
    for(int i = 0; i < num_coef; ++i)
    {
        assert(coe_res[i] == coef[i]);
    }
    return 0;
}