from nclib import Netcat
from Crypto.Util.number import long_to_bytes

#####################
# Functions to get  #
# stuff from server #
#####################

def ru(nc, b):
    """recv until"""
    r = b""
    while b not in r:
        r += nc.recv(1)
    return r

def get_pub_key():
    """Get public key"""
    ru(nc, b"n: ")
    n = int(nc.recvline())
    e = 3
    return n,e

def get_enc(buf):
    """get encryption of `buf` from server"""
    ru(nc, b"opt: ")
    nc.send(b"1\n")
    ru(nc, b"msg: ")
    nc.send(str(buf).encode() + b"\n")
    ru(nc, b"c: ")
    return int(nc.recvline())

def get_encflag():
    """get encryption of flag from server"""
    ru(nc, b"opt: ")
    nc.send(b"2\n")
    ru(nc, b"c: ")
    return int(nc.recvline())

def gen_poly(pad, padbitlength, ct, pt):
    """
    Generate polynomial enc(pad(pt)) - ct
    """
    return ((2<<(8*254)) + (pad * 2^(8*(256 - padbitlength//8 - 2))) + pt)^3 - ct

nc = Netcat(("193.57.159.27", 56926))
n,e = get_pub_key()

pt = 1 << (1920) # Long enough pt such that padding is 12 bytes

pads = []
for i in range(624//3):
    ct = get_enc(pt)
    x = PolynomialRing(Zmod(n), 'x').gen() # `x` represents padding
    pol = gen_poly(x, 96, ct, pt).monic() # create polynomial in `x`
    pads.append(pol.small_roots()[0]) # solve polynomial to get `x`
    
enc_flag = get_encflag() # Get encryption of the flag

########################
# Cloning the RNG used #
# for the padding and  #
# recovering the flag  #
########################

# https://github.com/JuliaPoo/MT19937-Symbolic-Execution-and-Solver
from MT19937 import MT19937

data624 = [] # 624 outputs of the server's MT19937
for p in pads:
    for i in range(3):
        data624.append((int(p) >> 32*i) & ((1<<32) - 1))

def copy_genrandbits(rng, nbits):
    """getrandbits implementation"""
    ds = [rng() for _ in range(nbits//32)][::-1]
    res = ds[0]
    for d in ds[1:]:
        res <<= 32
        res += d
    q = nbits % 32
    if q:
        res += (rng() >> (32-q)) << (32*(nbits//32))
    return res

for flen in range(54,0,-1): # Gues flag length (flen)
    
    padbitlength = (256-flen-3)*8

    rng_clone = MT19937(state_from_data = (data624, 32))
    for i in range(624):
        rng_clone()
    fpad = copy_genrandbits(rng_clone, padbitlength) # Predicted padding
    
    x = PolynomialRing(Zmod(n), 'flag').gen() # `x` represents the flag
    pol = gen_poly(fpad, padbitlength, enc_flag, x)
    roots = pol.small_roots(X=(1<<(flen*8))) # Solve for `x`

    if len(roots) != 0: # Solution found!
        break
        
flag = long_to_bytes(roots[0])
print("Flag:", flag.decode())