import time
import numpy as np
t = time.time()

import mini as orig # Import challenge code
from mini import *  # Import challenge code
from out import *   # Import output of challenge code given to us

F = GF(3) # Field we are working in

######################
# Process data given #
######################

# Convert known ct to the 3x3 tryte state
server_ct = [up(int_to_tyt(byt_to_int(ct)), W_SIZE ** 2, int_to_tyt(0)[0])[-1] for ct in all_enc]
# Flatten the ct to 27 trit vectors
server_ctv = [vector(F, [i for j in ct for i in j]) for ct in server_ct]
# Convert known pt to the 3x3 tryte state
server_pt = [up(int_to_tyt(pt), W_SIZE ** 2, int_to_tyt(0)[0])[-1] for pt in range(3^9)]
# Flatten the pt to 27 trit vectors
server_ptv = [vector(F, [i for j in pt for i in j]) for pt in server_pt]


#############################
# Create symbolic variables #
#############################

ptkr = PolynomialRing(F, ['pt%d'%i for i in range(9*3)] + ['kr%d'%i for i in range(9*3*8)]).gens()

# Flattened pt variables
pt_v = ptkr[:9*3]
# pt variables in 3x3 trit format
pt = tuple(tuple(pt_v[i*3:i*3+3]) for i in range(9))
# Flattened kr variables
kr_v = ptkr[9*3:]
# kr variables in 3x3 trit round keys format
kr = tuple(tuple(tuple(kr_v[i*9*3:i*9*3+9*3][j*3:j*3+3]) for j in range(9)) for i in range(8))

print("[*] Created symbolic variables!")

###################################
# Modify A3S to be GF(3) friendly #
###################################

xor = lambda a,b: a+b
uxor = lambda a,b: a-b
t_xor = lambda a,b: tuple(x+y for x,y in zip(a,b))
T_xor = lambda a,b: tuple(t_xor(i,j) for i,j in zip(a,b))
t_uxor = lambda a,b: tuple(x-y for x,y in zip(a,b))
T_uxor = lambda a,b: tuple(t_uxor(i,j) for i,j in zip(a,b))

SBOX_TYT = dict((int_to_tyt(i)[0], int_to_tyt(s)[0]) for i,s in enumerate(SBOX))
ISBOX = tuple(SBOX.index(i) for i in range(27))
ISBOX_TYT = dict((int_to_tyt(i)[0], int_to_tyt(s)[0]) for i,s in enumerate(ISBOX))

def sbox_affine(i:tuple):
    return (
        (0 + i[0]*1 + i[1]*1 + i[2]*0),
        (0 + i[0]*0 + i[1]*0 + i[2]*1),
        (1 + i[0]*0 + i[1]*2 + i[2]*2)
    )

def expand(tyt):
    words   = tyt_to_wrd(tyt) 
    size    = len(words)
    rnum    = size + 3
    rcons   = rcon(rnum * 3 // size)

    for i in range(size, rnum * 3):
        k   = words[i - size]
        l   = words[i - 1]
        if i % size == 0:
            s = tuple(sbox_affine(i) for i in rot_wrd(l))
            k = T_xor(k, s)
            k = (t_xor(k[0], rcons[i // size - 1]),) + k[1:]
        else:
            k = T_xor(k, l)
        words = words + (k,)

    return up(down(words[:rnum * 3]), W_SIZE ** 2, int_to_tyt(0)[0])

def tri_mulmod(A, B, mod=POLY):
    c = [0] * (len(mod) - 1)
    for a in A[::-1]:
        c = [0] + c
        x = tuple(b * a for b in B)
        c[:len(x)] = t_xor(c, x)
        n = -c[-1]*mod[-1]
        c[:] = [x+y*n for x,y in zip(c,mod)]
        c.pop()
    return tuple(c)

def tyt_mulmod(A, B, mod=POLY2, mod2=POLY):
    fil = [(0,) * T_SIZE]
    C = fil * (len(mod) - 1)
    for a in A[::-1]:
        C = fil + C
        x = tuple(tri_mulmod(b, a, mod2) for b in B)
        C[:len(x)] = T_xor(C, x)
        num = modinv(mod[-1], mod2)
        num2 = tri_mulmod(num, C[-1], mod2)
        x = tuple(tri_mulmod(m, num2, mod2) for m in mod)
        C[:len(x)] = T_uxor(C, x)

        C.pop()
    return C

def add(a,b):
    return tuple(
        tuple(x+y for x,y in zip(i,j)) for i,j in zip(a,b)
    )

def sub(a):
    return tuple(
        sbox_affine(x) for x in a
    )

def shift(a):
    return [
        a[i] for i in SHIFT_ROWS
    ]

def mix(tyt):
    tyt = list(tyt)
    for i in range(W_SIZE):
        tyt[i::W_SIZE] = tyt_mulmod(tyt[i::W_SIZE], CONS)
    return tuple(tyt)

    

#########################################
# Create constraints from key expansion #
#########################################

# xkey_mat1 * kr = xkey_const1
# xkey_mat2 * kr = xkey_const2 (except for maybe ~3 constraints where sub isn't affine)

xkey_mat1 = []   # An array of 135 x 216 numbers
xkey_const1 = [] # A vector of 135 numbers

for i in range(6,24):
    if i%5 == 0: continue
    for j in range(9):
        r = vector([0]*3*8*9)
        r[i*9+j] = 1; r[i*9-9+j] = 2; r[i*9-9*5+j] = 2
        xkey_mat1.append(r)
        xkey_const1.append(0)

# Symbolically execute the 2nd line of the key expansion
# and save the eqns
xkey_eqns = []
for i in range(5,25,5):
    rcons = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 2))
    k  = [tuple(kr_v[(i-5)*9:(i-5)*9+9][j*3:j*3+3]) for j in range(3)]
    l  = [tuple(kr_v[(i-1)*9:(i-1)*9+9][j*3:j*3+3]) for j in range(3)]
    s = sub(rot_wrd(l))
    k0 = T_xor(k, s)
    k0 = (t_xor(k0[0], rcons[i // 5 - 1]),) + k0[1:]
    k0 = [i for j in k0 for i in j]
    k0 = [i-j for i,j in zip(k0,kr_v[i*9:i*9+9])]
    xkey_eqns.extend(k0)

# Create matrix from the eqns above
xkey_mat2 = []
xkey_const2 = []
for k in xkey_eqns:
    r = vector([0]*3*8*9)
    s = 0
    for v,c in k.dict().items():
        if 1 not in v:
            s -= c
            continue
        else:
            vi = list(v).index(1)
            r[vi - 9*3] += c
    xkey_mat2.append(r)
    xkey_const2.append(s)

print("[*] Created constraints from key expansion!")

#######################################
# Create constraints from pt-ct pairs #
#######################################

def forward_to_checkpoint(ctt, keys, checkpoint):

    """
    Encrypts ctt with keys to the checkpoint
        ctt is the plaintext in 3x3 tryte format
        keys is the expanded key
    """

    ctt = add(ctt, keys[0])

    for r in range(1, len(keys) - 1):
        ctt = sub(ctt)
        ctt = shift(ctt)
        ctt = mix(ctt)
        ctt = add(ctt, keys[r])

    if checkpoint==0: return ctt
    ctt  = sub(ctt)
    ctt  = shift(ctt)
    if checkpoint==1: return ctt
    ctt  = add(ctt, keys[-1])

    return ctt

def gen_mat(checkpoint):

    """
    Generates matrices and vectors corresponding
    to the affine transform to `checkpoint`.

    Returns ptmat, ptconst, kptmat where
        ct = ptmat*pt + ptconst + kptmat*kr
    """
    
    # pt and kr are symbolic
    ctt = forward_to_checkpoint(pt, kr, checkpoint)

    # Create the matrix from resulting
    # symbolic equation
    ptmat = []
    kptmat = []
    ptconst = []
    for w in ctt:
        for t in w:
            rpt = vector([0]*9*3)
            rkpt = vector([0]*9*3*8)
            s = 0
            for v,c in t.dict().items():
                if 1 not in v: 
                    s += c; continue
                vi = list(v).index(1)
                if vi>=9*3: 
                    rkpt[vi-9*3] += c
                else:    
                    rpt[vi] += c
            ptmat.append(rpt)
            ptconst.append(s)
            kptmat.append(rkpt)
    
    # ct = ptmat*pt + ptconst + kptmat*kr
    return matrix(F, ptmat), vector(F, ptconst), matrix(F, kptmat)


# -------------#
# Checkpoint 1 #
# -------------#

print("[*] Exploiting checkpoint 1...")

# Calculate `ptmat`
ptmat, _, _ = gen_mat(0)
# Encrypt all known pt up till checkpoint 0
spt = ptmat * matrix(F, server_ptv).T
spt = spt.T

# Collect 3^9*8 pt differences
# Offsets used are just powers of 3
# No reason for it, it just looks nice
dspt = []
for offset in [1,3,9,27,81,243,729,2187]:
    dspt += [(spt[i]-spt[(i+offset)%3^9], (i,(i+offset)%3^9)) for i in range(3^9)]

# An array of 9 dictionaries
# > each element corresponds to each tryte of the last_key
# > each element is a dictionary containing the "score" of each
#   possible guess
# > The higher the score the more probable the tryte is the 
#   correct guess
all_kscore = []

for cidx in range(9): # enumerate all trytes of `last_key`

    pidx = SHIFT_ROWS[cidx]
    kscore = [0]*27
    for dptidx in range(len(dspt)):

        dpt,(i0,i1) = dspt[dptidx]
        dct = (server_ct[i0], server_ct[i1])

        c = (dct[0][cidx], dct[1][cidx])
        p = tuple(dpt[pidx*3:pidx*3+3])

        for k in range(27): # enumerate all tryte guesses
 
            # Partial decrypt
            kt = orig.int_to_tyt(k)[0]
            ci = (orig.t_uxor(c[0],kt), orig.t_uxor(c[1],kt)) # unadd
            ci = (ISBOX_TYT[ci[0]], ISBOX_TYT[ci[1]]) # unsub
            ci = orig.t_uxor(ci[0], ci[1])

            if ci == p: # if matches, add 1 to score
                kscore[k] += 1

        print(dptidx, end="\r")
        
    all_kscore.append(kscore)
    print("[-]", cidx, 'done!')

# Get the k with the highest scores as the last_key
last_key = [int_to_tyt(all_kscore[i].index(max(all_kscore[i])))[0] for i in range(9)]
last_key = [i for j in last_key for i in j]

print("[*] Exploited checkpoint 1!")
print("[#] Last round key:", last_key)


# -------------#
# Checkpoint 2 #
# -------------#

print("[*] Exploiting checkpoint 2...")

ptmat, ptconst, kptmat = gen_mat(1)

# Get all 3^9 values of ct-A(pt)
spt = ptmat * matrix(F, server_ptv).T + matrix(F, [list(ptconst)]*3^9).T
rhs_prob = matrix(F, server_ctv).T - spt

# Get most probable trit
# rhs = ct - A(pt)
from collections import Counter
rhs = [sorted(Counter(i).items(), key=lambda x:x[1])[-1][0] for i in rhs_prob]

# Calculate M
# M*kr = rhs should hold
M = np.zeros((27,216))
for i in range(27):
    M[i][216-27+i] = 1
M = matrix(F,M) + kptmat

print("[*] Exploited checkpoint 2!")


###########################
# Putting it all together #
###########################

print("[*] Putting everything together!")

def build_solve_mat(xkey_mat2, xkey_const2):

    """
    Collate all linear constraints together.
    > arguments corresponds to the constraints from the
      key expansion that requires SBOX. We might need to
      remove some of those constraints, hence they are placed
      as arguments for easy modifications
    """
    
    # Form the 54 constraints
    
    # 27 from knowing the last_key
    a3s_lhs = [[0]*216 for _ in range(27)]
    a3s_rhs = last_key.copy()
    for i in range(27):
        a3s_lhs[i][216-27+i] = 1
        
    # 27 from assuming the whole cipher is linear
    a3s_lhs.extend(list(M))
    a3s_rhs.extend(rhs)

    # Combine with key expansion
    
    # From expansion that doesn't involve sbox
    a3s_lhs.extend(list(xkey_mat1))
    a3s_rhs.extend(list(xkey_const1))
    
    # From expansion that involves sbox (that might be wrong)
    a3s_lhs.extend(list(xkey_mat2))
    a3s_rhs.extend(list(xkey_const2))

    a3s_lhs = matrix(F, a3s_lhs)
    a3s_rhs = vector(F, a3s_rhs)
    return a3s_lhs, a3s_rhs


# Assume one sub wasnt affine in the key expansion
for i in range(12): # Try all 12 possible substitutions
    a = xkey_mat2.copy()
    b = xkey_const2.copy()
    a = [p for j,p in enumerate(a) if j not in range(i*3,i*3+3)]
    b = [p for j,p in enumerate(b) if j not in range(i*3,i*3+3)]
    l,r = build_solve_mat(a,b)
    try:
        key_recovered = l.solve_right(r)
        print("[*] kr found! ^-^")
        break
    except: # Crashes if l.solve_right has no solutions
        pass
    
# Check if the equation l*kr = r has more than 1 solution
assert len(l.right_kernel().basis()) == 0, "Not the only solution!"

key_recovered = [int(i) for i in key_recovered]
key_recovered = tyt_to_int(tuple(tuple(key_recovered[i*3:i*3+3]) for i in range(15)))
key_recovered = int_to_byt(key_recovered)
print("[*] Original Key:", key_recovered.hex())

from hashlib import sha512

hsh = sha512(key_recovered).digest()
flag = byte_xor(hsh, enc_flag)
print("[*] Flag:", flag.decode())
print("[*] Time taken: %ds"%(time.time()-t))