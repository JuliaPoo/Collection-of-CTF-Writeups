from a3s import * # import the challenge code

##########################################################
# Modify some of the functions to be sage GF(3) friendly #
##########################################################

t_xor = lambda a,b: tuple(x+y for x,y in zip(a,b))
T_xor = lambda a,b: tuple(t_xor(i,j) for i,j in zip(a,b))
t_uxor = lambda a,b: tuple(x-y for x,y in zip(a,b))
T_uxor = lambda a,b: tuple(t_uxor(i,j) for i,j in zip(a,b))

def sbox_affine(i:tuple):
    return (
        (0 + i[0]*1 + i[1]*2 + i[2]*1),
        (2 + i[0]*0 + i[1]*1 + i[2]*2),
        (0 + i[0]*2 + i[1]*1 + i[2]*0)
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

#####################
# Symbolise the key #
#####################

F = GF(3)
FP = PolynomialRing(F, 'k', 25*3*3)
keyFP = FP.gens()

##################
# Expand the key #
##################

exkeyFP = tuple(tuple(keyFP[i*3+j] for j in range(3)) for i in range(25*3))
exkeyFP = expand(exkeyFP)

################################
# Perform encryption of `sus.` #
# with symbolised key          #
################################

m = byt_to_int(b'sus.')
m = up(int_to_tyt(m), W_SIZE ** 2, int_to_tyt(0)[0])[-1]

assert len(exkeyFP) == KEYLEN

ctt = add(m, exkeyFP[0])

for r in range(1, len(exkeyFP) - 1):
    ctt = sub(ctt)
    ctt = shift(ctt)
    ctt = mix(ctt)
    ctt = add(ctt, exkeyFP[r])

ctt  = sub(ctt)
ctt  = shift(ctt)
ctt  = add(ctt, exkeyFP[-1])

###############################
# Form linear equation with   #
# enc(b"sus.") as RHS         #
###############################

mat = []
const = []
for x in ctt:
    for y in x:
        s = vector([0]*75*3)
        c = 0
        for i,j in y.dict().items():
            if sum(i) == 0:
                c += j
            s += vector(i)*j
        const.append(c)
        mat.append(s)
mat = Matrix(F, mat)

rhs = int_to_tyt(byt_to_int(b'\x06\x0f"\x02\x8e\xd1'))
rhs = vector(F, [i for j in rhs for i in j])
rhs -= vector(F, const)

#############################
# Solve the linear equation #
# to get key                #
#############################

key_rec = mat.solve_right(rhs)
key_rec = tuple(tuple(key_rec[3*i+j] for j in range(3)) for i in range(75))

print(key_rec)

