#include "linear_gkr/prime_field.h"
#include <climits>
#include <ctime>
namespace prime_field
{

	const __uint128_t __max_ull = ULLONG_MAX;
	const __uint128_t __zero = 0LL;
	const __uint128_t __max_128 = ~(__zero);



	
	u256b::u256b(){lo = mid = 0; hi = 0;}
	u256b::u256b(const unsigned long long &x)
	{
		lo = x;
		mid = 0;
		hi = 0;
	}
	u256b::u256b(const __uint128_t &x)
	{
		lo = (unsigned long long)(x & __max_ull);
		mid = (unsigned long long)(x >> 64);
		hi = 0;
	}

	u256b::u256b(const char* x, int len, int base)
	{
		u256b ret = (u256b)(unsigned long long)0;
		u256b tmp = (u256b)(unsigned long long)1;
		for(int i = len - 1; i >= 0; --i)
		{
			ret = ret + tmp * ((u256b)(unsigned long long)(x[i] - '0'));
			tmp = tmp * ((u256b)(unsigned long long)base);
		}
		*this = ret;
	}

	inline u256b u256b::operator + (const u256b &x) const
	{
		u256b ret;
		bool carry;
		__uint128_t midd;
		ret.lo = lo + x.lo;
		carry = (ret.lo < lo);
		midd = (__uint128_t)mid + (__uint128_t)x.mid;
		if(midd == __max_ull)
			midd += carry;
		ret.mid = mid + x.mid + carry;
		ret.hi = hi + x.hi + (midd >> 64);
		return ret;
	}

	inline u256b u256b::operator - (const u256b &x) const
	{
		u256b not_x;
		not_x.lo = ~x.lo;
		not_x.mid = ~x.mid;
		not_x.hi = ~x.hi;
		return *this + not_x;
	}

	u256b u256b::operator * (const u256b &x) const
	{
		u256b ret;
		__uint128_t lolo = (__uint128_t)lo * (__uint128_t)x.lo;
		
		__uint128_t lomid1 = (__uint128_t)mid * (__uint128_t)x.lo;
		__uint128_t lomid2 = (__uint128_t)lo * (__uint128_t)x.mid;
		
		//hi * hi omitted
		ret.lo = (unsigned long long)lolo;
		__uint128_t carry = lolo >> 64; //this carry is less than 2**64
		ret.mid = ((unsigned long long)lomid1 + (unsigned long long)lomid2 + (unsigned long long)carry);
		//this carry is not necessary less than 2**64
		__uint128_t carry2 = ((lomid1 >> 64) + (lomid2 >> 64)) + 
						   (((lomid1 & __max_ull) + (lomid2 & __max_ull) + carry) >> 64);
						   
		ret.hi = (__uint128_t)lo * (__uint128_t)x.hi + (__uint128_t)hi * (__uint128_t)x.lo + 
				 (__uint128_t)mid * (__uint128_t)x.mid + (((__uint128_t)mid * (__uint128_t)x.hi + (__uint128_t)hi * (__uint128_t)x.mid) << 64) + carry2;
		return ret;
	}
	/*
	inline u256b u256b::operator << (const int &x) const{
		u256b ret;
		for(int i = 0; i < x; i++)
			ret << 1;
		ret = ret % prime_field::mod;
		return ret;
	}

	inline u256b u256b::operator >> (const int &x) const{
		u256b ret;
		for(int i = 0; i < x; i++)
			ret >> 1;
		ret = ret % prime_field::mod;
		return ret;
	}
	*/
	inline u256b u256b::left_128()
	{
		u256b ret;
		ret.lo = 0;
		ret.mid = 0;
		ret.hi = ((__uint128_t)lo | ((__uint128_t)mid << 64));
		return ret;
	}
	
	inline u256b u256b::operator & (const u256b &x) const
	{
		u256b ret;
		ret.lo = lo & x.lo;
		ret.mid = mid & x.mid;
		ret.hi = hi & x.hi;
		return ret;
	}
	inline int u256b::bit_at(int pos) const
	{
		if(pos < 64)
			return (lo >> pos) & 1;
		if(pos < 128)
			return (mid >> (pos - 64)) & 1;
		else
			return (hi >> (pos - 128)) & 1;
	}
	inline bool u256b::operator <= (const u256b &x) const
	{
		if(hi < x.hi)
			return true;
		if(hi > x.hi)
			return false;
		if(mid < x.mid)
			return true;
		if(mid > x.mid)
			return false;
		if(lo <= x.lo)
			return true;
		return false;
	}
	inline bool u256b::operator != (const u256b &x) const
	{
		return hi != x.hi || lo != x.lo || mid != x.mid;
	}

	inline bool u256b::operator > (const u256b &x) const
	{
		if(hi > x.hi)
			return true;
		if(hi < x.hi)
			return false;
		if(mid > x.mid)
			return true;
		if(mid < x.mid)
			return false;
		return lo > x.lo;
	}
	inline bool u256b::operator < (const u256b &x) const
	{
		if(hi < x.hi)
			return true;
		if(hi > x.hi)
			return false;
		if(mid < x.mid)
			return true;
		if(mid > x.mid)
			return false;
		return lo < x.lo;
	}

	inline u256b midhi_mul(const u256b &a, const u256b &b)
	{
		u256b ret;
		__uint128_t lolo = (__uint128_t)a.lo * (__uint128_t)b.lo;
		
		__uint128_t lomid1 = (__uint128_t)a.mid * (__uint128_t)b.lo;
		__uint128_t lomid2 = (__uint128_t)a.lo * (__uint128_t)b.mid;
		
		//hi * hi omitted
		ret.lo = (unsigned long long)lolo;
		__uint128_t carry = lolo >> 64; //this carry is less than 2**64
		ret.mid = ((unsigned long long)lomid1 + (unsigned long long)lomid2 + (unsigned long long)carry);				   
		return ret;
	}

	u256b interm_tmp;
	u256b interesting_tmps[256];
	u256b interesting_combination[256 / 8][1 << 8];
	u256b one, zero;
	void preprocess_mod(u256b mod)
	{
		u256b not_mod;
		not_mod.lo = ~mod.lo;
		not_mod.mid = ~mod.mid;
		not_mod.hi = ~mod.hi;
		not_mod = not_mod + (u256b)(unsigned long long)1;
		interesting_tmps[0] = interm_tmp;
		for(int i = 1; i < 256; ++i)
		{
			interesting_tmps[i] = interesting_tmps[i - 1] + interesting_tmps[i - 1];
			if(mod <= interesting_tmps[i])
			{
				interesting_tmps[i] = interesting_tmps[i] + not_mod;
			}
		}

		for(int i = 0; i < 256 / 8; ++i)
		{
			for(int j = 0; j < (1 << 8); ++j)
			{
				u256b sum = (unsigned long long)0;
				for(int k = 0; k < 8; ++k)
				{
					if((j >> k) & 1)
					{
						sum = sum + interesting_tmps[i * 8 + k];
						if(mod <= sum)
							sum = sum + not_mod;
					}
				}
				interesting_combination[i][j] = sum;
			}
		}
	}

	u512b::u512b(const u256b &x)
	{
		lo = ((__uint128_t)x.mid << 64) | (__uint128_t)x.lo;
		mid = x.hi;
		hi.lo = 0;
		hi.mid = 0;
		hi.hi = 0;
	}
	u512b::u512b(const __uint128_t &x)
	{
		lo = x;
		mid = 0;
		hi.lo = 0;
		hi.mid = 0;
		hi.hi = 0;
	}
	u512b::u512b(const char* x, int len, int base)
	{
		u512b ret = (u256b)(unsigned long long)0;
		u512b tmp = (u256b)(unsigned long long)1;
		for(int i = 0; i < len; ++i)
		{
			ret = ret + tmp * (u512b)((u256b)(unsigned long long)(x[i] - '0'));
			tmp = tmp * (u512b)((u256b)(unsigned long long)10);
		}
		*this = ret;
	}

	u512b::u512b(){lo = mid = 0; hi.lo = 0; hi.mid = 0; hi.hi = 0;}
	u512b u512b::operator + (const u512b &x) const
	{
		u512b ret;
		__uint128_t carry, carry2;
		ret.lo = lo + x.lo;
		carry = ret.lo < lo;
		ret.mid = mid + x.mid + carry;
		if(carry == 0)
			carry2 = ret.mid < mid;
		else
			carry2 = ret.mid <= mid;
		ret.hi = hi + x.hi + carry2;
		return ret;
	}

	u512b u512b::operator - (const u512b &x) const
	{
		u512b not_x;
		not_x.hi.hi = ~x.hi.hi;
		not_x.hi.mid = ~x.hi.mid;
		not_x.hi.lo = ~x.hi.lo;
		not_x.mid = ~x.mid;
		not_x.lo = ~x.lo;
		not_x = not_x + one;
		return *this + not_x;
	}
	
	u512b u512b::operator * (const u512b &x) const
	{
		u512b ret;
		u256b lolo = (u256b)lo * (u256b)x.lo;
		
		u256b lomid1 = (u256b)mid * (u256b)x.lo;
		u256b lomid2 = (u256b)lo * (u256b)x.mid;
		
		u256b lohi1 = (u256b)lo * (u256b)x.hi;
		u256b lohi2 = (u256b)hi * (u256b)x.lo;
		u256b midmid = (u256b)mid * (u256b)x.mid;
		
		//u256b midhi1 = (u256b)mid * (u256b)x.hi;
		//u256b midhi2 = (u256b)hi * (u256b)x.mid;
		
		u256b midhi1 = midhi_mul(mid, x.hi);
		u256b midhi2 = midhi_mul(hi, x.mid);
		
		//hi * hi omitted
		ret.lo = (__uint128_t)lolo.lo | ((__uint128_t)lolo.mid << 64);
		__uint128_t carry = lolo.hi; //this carry is less than 2**128
		
		u256b tmp = (lomid1 + lomid2 + carry);
		
		ret.mid = (__uint128_t)tmp.lo | ((__uint128_t)tmp.mid << 64);
		//this carry is not necessary less than 2**128
		u256b carry2 = ((u256b)(lomid1.hi) + (u256b)(lomid2.hi)) + 
						   (((u256b)((__uint128_t)lomid1.lo | ((__uint128_t)lomid1.mid << 64))
						   + (u256b)((__uint128_t)lomid2.lo | ((__uint128_t)lomid2.mid << 64)) + (u256b)carry).hi);
						   
		ret.hi = lohi1 + lohi2 + midmid + ((midhi1 + midhi2).left_128()) + carry2;
		return ret;
	}
	
	u512b u512b::operator % (const u256b &x) const
	{
		u256b ret = zero;
		u256b not_x;
		not_x.lo = ~x.lo;
		not_x.mid = ~x.mid;
		not_x.hi = ~x.hi;
		not_x = not_x + one;
		u256b lowbits = (u256b)lo + ((u256b)mid).left_128();
		while(x <= lowbits)
			lowbits = lowbits + not_x;
		ret = lowbits;
		int index[32];
		for(int i = 0; i < 32; ++i)
		{
			if(i < 8)
			{
				index[i] = (hi.lo >> (i * 8)) & ((1 << 8) - 1);
			}
			else if(i < 16)
			{
				index[i] = (hi.mid >> (i * 8 - 64)) & ((1 << 8) - 1);
			}
			else
			{
				index[i] = (hi.hi >> (i * 8 - 128)) & ((1 << 8) - 1);
			}
		}
		for(int i = 0; i < 8; ++i)
		{
			u256b s0, s1;
			s0 = interesting_combination[i * 4][index[i * 4]] + 
				 interesting_combination[i * 4 + 1][index[i * 4 + 1]];
			s1 = interesting_combination[i * 4 + 2][index[i * 4 + 2]] + 
				 interesting_combination[i * 4 + 3][index[i * 4 + 3]];
			if(x <= s0)
				s0 = s0 + not_x;
			if(x <= s1)
				s1 = s1 + not_x;
			u256b sum = s0 + s1;
			if(x <= sum)
				sum = sum + not_x;
			ret = ret + sum;
			if(x <= ret)
				ret = ret + not_x;
		}
		return u512b(ret);
	}
	bool u512b::operator != (const u512b &x) const
	{
		return lo != x.lo || mid != x.mid || hi != x.hi;
	}
	bool u512b::operator > (const u512b &x) const
	{
		if(hi > x.hi)
			return true;
		if(hi < x.hi)
			return false;
		if(mid > x.mid)
			return true;
		if(mid < x.mid)
			return false;
		return lo > x.lo;
	}

	bool u512b::operator >= (const u512b &x) const
	{
		if(hi > x.hi)
			return true;
		if(hi < x.hi)
			return false;
		if(mid > x.mid)
			return true;
		if(mid < x.mid)
			return false;
		return lo >= x.lo;
	}

	bool u512b::operator < (const u512b &x) const
	{
		if(hi < x.hi)
			return true;
		if(hi > x.hi)
			return false;
		if(mid < x.mid)
			return true;
		if(mid > x.mid)
			return false;
		return lo < x.lo;
	}
	void u512b::random()
	{
		lo = rand();
		mid = rand();
		hi = (u256b)(unsigned long long)rand();
	}



	u256b mod;
	bool initialized = false;
	//independent_bits_engine<mt19937, 256, cpp_int> gen;

	void init(std::string s, int base)
	{
		assert(base == 10);
		initialized = true;
		mod = u256b(s.c_str(), s.length(), base);
		std::string interm_tmp_s = "15003436851221201713926160155297504394712507045212047545287310822919708344242";
		interm_tmp = u256b(interm_tmp_s.c_str(), interm_tmp_s.length(), 10);
		preprocess_mod(mod);
		one = (u256b)(unsigned long long)1;
		zero = (u256b)(unsigned long long)0;
	}
	void init_random()
	{
	}
	field_element::field_element()
	{
		value = (u256b)(unsigned long long)0;
	}
	field_element::field_element(const int x)
	{
		value = (u256b)(unsigned long long)x;
	}
	field_element field_element::operator + (const field_element &b) const
	{
		field_element ret;
		ret.value = (b.value + value);
		if(ret.value >= mod)
			ret.value = ret.value - mod;
		return ret;
	}
	field_element field_element::mul_non_mod(const field_element &b) const
	{
		field_element ret;
		ret.value = (b.value * value);
		return ret;
	}
	field_element field_element::operator * (const field_element &b) const
	{
		field_element ret;
		ret.value = (b.value * value) % mod;
		return ret;
	}
	field_element field_element::operator / (const field_element &b) const
	{
		//todo
		assert(false);
		//return ret;
	}
	field_element field_element::operator - (const field_element &b) const
	{
		field_element ret;
		if(value >= b.value)
			ret.value = value - b.value;
		else
			ret.value = value + mod - b.value;
		return ret;
	}
	char* field_element::to_string()
	{
		static char ret[512];
		for(int i = 0; i < 128; ++i)
			ret[i] = ((value.lo >> i) & 1) + '0';
		for(int i = 0; i < 128; ++i)
			ret[i + 128] = ((value.mid >> i) & 1) + '0';
		for(int i = 0; i < 256; ++i)
			ret[i + 256] = ((value.hi.bit_at(i)) & 1) + '0';
		return ret;
	}
	field_element random()
	{
		field_element ret;
		ret.value.random();
		ret.value = ret.value % mod;
		return ret;
	}
	bool field_element::operator != (const field_element &b) const
	{
		return value != b.value;
	}
}
