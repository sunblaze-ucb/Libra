#ifndef __transcript
#define __transcript
#include "linear_gkr/prime_field.h"
#include <vector>
#include "bn.h"
#include "infrastructure/my_hhash.h"
class container
{
public:
    static const int fld_ele_indicator = 1;
    static const int curve_ele_indicator = 0;
    static const int randomness_indicator = 2;
    int indicator;
    prime_field::field_element fld_ele;
    bn::Ec1 curve_ele;
    container(){}
    container(prime_field::field_element x, int ind)
    {
        fld_ele = x;
        indicator = ind;
    }
    container(bn::Ec1 x)
    {
        curve_ele = x;
        indicator = curve_ele_indicator;
    }
    std::vector<bool> output()
    {
        //compress and output
        std::vector<bool> out;
        if(indicator == fld_ele_indicator || indicator == randomness_indicator)
        {
            auto x = fld_ele.bit_stream();
            out.insert(out.end(), x.begin(), x.end());
        }
        else
        {
            const int bytes = sizeof(bn::Ec1);
            char* x = (char*)&curve_ele;
            for(int i = 0; i < bytes; ++i)
            {
                for(int j = 0; j < 8; ++j)
                    out.push_back((x[i] >> j) & 1);
            }
        }
        return out;
    }
};
class proof_transcript
{
public:
    std::vector<container> msg;
    std::vector<bool> output_bit_stream(bool need_rand)
    {
        std::vector<bool> output;
        for(auto x : msg)
        {
            if(x.indicator == x.randomness_indicator && !need_rand)
                continue;
            auto out = x.output();
            output.insert(output.end(), out.begin(), out.end());
        }
        return output;
    }
    void output_to_file(char* file)
    {
        FILE *f = fopen(file, "wb");
        auto vec = output_bit_stream(false);
        int bits = vec.size();
        int bytes;
        if(bits % 8 == 0)
            bytes = bits / 8;
        else
            bytes = bits / 8 + 1;
        
        char* buf = new char[bytes];
        for(int i = 0; i < bytes; ++i)
            buf[i] = 0;
        for(int i = 0; i < bits; ++i)
        {
            buf[i / 8] |= ((char)vec[i] << (i & 7));
        }
        fwrite(buf, sizeof(char), bytes, f);
        delete[] buf;
        fclose(f);
    }

    prime_field::field_element* random(int size)
    {
        auto x = output_bit_stream(true);
        //align with 512
        int data[512 / (8 * sizeof(int))];
        while(x.size() % 512 != 0)
            x.push_back(0);
        for(int i = 0; i < (512 / (8 * sizeof(int))); ++i)
            data[i] = 0;
        for(int i = 0; i < 512; ++i)
            data[i / 32] |= (((int)x[i]) << (i & (8 * sizeof(int) - 1)));
        __hhash_digest h;
        my_hhash(data, &h);
        int cur = 512;
        while(cur < x.size())
        {
            for(int i = 0; i < (512 / (8 * sizeof(int))); ++i)
                data[i] = 0;
            memcpy(data, &h, sizeof(h));
            for(int i = 0; i < 256; ++i)
            {
                int curbit = x[cur++];
                data[(i + 256) / 32] |= curbit << ((i + 256) & 31);
            }
            my_hhash(data, &h);
        }
        prime_field::field_element* r;
        r = new prime_field::field_element[size];
        for(int i = 0; i < size; ++i)
        {
            memcpy(&r[i].value.lo, &h.h0, sizeof(r[i].value.lo));
            memcpy(&r[i].value.mid, &h.h1, sizeof(r[i].value.mid));
            r[i].value.hi.lo = r[i].value.hi.mid = 0;
            r[i].value.hi.hi = 0;
            r[i].value = r[i].value % prime_field::mod;
            for(int i = 0; i < (512 / (8 * sizeof(int))); ++i)
                data[i] = 0;
            memcpy(data, &h, sizeof(h));
            my_hhash(data, &h);
        }
        for(int i = 0; i < size; ++i)
        {
            msg.push_back(container(r[i], container::randomness_indicator));
        }
        return r;
    }
    prime_field::field_element random()
    {
        prime_field::field_element r;
        auto r_arr = random(1);
        r = r_arr[0];
        delete[] r_arr;
        return r;
    }
    int get_proof_size()
    {
        int ret = 0;
        for(int i = 0; i < msg.size(); ++i)
        {
            if(msg[i].indicator != container::randomness_indicator)
            {
                if(msg[i].indicator == container::fld_ele_indicator)
                    ret += sizeof(prime_field::field_element) / 2;
                else
                    ret += sizeof(bn::Ec1);
            }
        }
        return ret;
    }
};

#endif