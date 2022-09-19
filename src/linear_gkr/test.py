p = 16798108731015832284940804142231733909759579603404752749028378864165570215949
mod = p
shift = 254 * 2
factor = (1 << shift) // p

print shift
print factor
print factor.bit_length()

def reduce(x):
    t = x - (((x * factor) >> shift) * mod)
    if(t >= mod):
        return t - mod
    return t

#random test

r = 187927128565899168703584974214286520196 + (145534863299331107206045469267347456661 << 128) + (507605208255225040 << 256)


print reduce(r)
print r % mod