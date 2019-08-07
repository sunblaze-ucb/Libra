#ifndef __transcript
#define __transcript
#include "linear_gkr/prime_field.h"
#include <vector>
#include "bn.h"
class container
{
public:
    bool indicator;
    prime_field::field_element fld_ele;
    bn::Ec1 curve_ele;
    std::vector<bool> output()
    {
        const int fld_ele_indicator = 1;
        const int curve_ele_indicator = 0;
        //compress and output
        std::vector<bool> out;
        out.push_back(indicator);
        if(indicator == fld_ele_indicator)
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
    std::vector<bool> output_bit_stream()
    {
        std::vector<bool> output;
        for(auto x : msg)
        {
            auto out = x.output();
            output.insert(output.end(), out.begin(), out.end());
        }
        return output;
    }
    void output_to_file(char* file)
    {
        FILE *f = fopen(file, "wb");
        auto vec = output_bit_stream();
        int bits = vec.size();
        int bytes;
        if(bits % 8 == 0)
            bytes = bits / 8;
        else
            bytes = bits / 8 + 1;
        
        char* buf = new char[bytes];
        fwrite(buf, sizeof(char), bytes, f);
        delete[] buf;
        fclose(f);
    }

    prime_field::field_element random()
    {
        auto x = output_bit_stream();
        
    }
};

#endif