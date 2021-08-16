# The Challenge

The challenge is similar to the challenge [A3S](../A3S). We are given an implementation of a cipher `A3S` and a [.pdf](chal/help.pdf) file that briefly explains the cipher. A3S is very similar to AES, with some notable differences:

1. Trits are used instead of bits. This means working in `GF(3)` instead of `GF(2)`.
2. The key expansion allows for keys larger than a round key. This means recovering the last round key isn't enough to recover the original key. This will become relevant later.

We are also given `3^9` plaintext and ciphertext pairs, and we have to recover the original key used in order to get the flag. 

## Metadata

// TODO

# Some A3S Background

As mentioned before, A3S is very similar to [AES](https://en.wikipedia.org/wiki/Advanced_Encryption_Standard). A3S operates on a state of `27` trits, arranged as `3x3` trytes (`3` trits to a tryte, kinda like `8` bits to a byte). The plaintext is converted to this state, which then passes through some operations, and the final state is then converted to the ciphertext.

Just like AES, A3S has the operations for encryption `sub`, `shift_rows`, `mix_columns`, `add_round_key`. More details in the [challenge pdf](chal/help.pdf).

An important thing to note is that `shift_rows`, `mix_columns` and `add_round_key` are all linear operations on the state. In a secure implementation, the `sub` would be the crucial component ensuring the cipher isn't affine, which is Very Very Bad™. See my solution for the challenge [A3S](../A3S). The lookup table for this substitution is the all important `SBOX`.

The key expansion for A3S relies on `sub` as well to ensure the key expansion isn't affine. This will be important later.

# The Vuln

From now on, assume we are working with `GF(3)` unless stated otherwise. This means if you were to see say `2+2`, take modulo 3 and have `2+2=1`.

I shamefully didn't notice when solving [A3S](../A3S) but the author provided a hint in the [A3S source](chal/mini.py): A comment `Secure enough ig` right above the SBOX.

But yes! Just like it's child challenge, the SBOX is vulnerable. This can be seen by plotting it's difference distribution table, which is plotting `a+b` over `SBOX[a] + SBOX[b]`.

```python
from mini import * # Import challenge code

import matplotlib.pyplot as plt
import numpy as np

mat = np.zeros((27,27))
for a in range(27):
    for b in range(27):
        sa, sb = SBOX[a], SBOX[b]
        i = tri_to_int(t_xor(int_to_tyt(a )[0], int_to_tyt(b )[0]))
        o = tri_to_int(t_xor(int_to_tyt(sa)[0], int_to_tyt(sb)[0]))
        mat[o][i] += 1
        
plt.imshow(mat)
plt.show()
```

<p align="center">
  <img  src="rsrc/sbox_diff.png" alt="SBOX differential table">
</p>

This looks bad, as a good SBOX should have a pretty uniform differential table. This is a good indication that the SBOX is _really_ close to being affine. It's pretty straightforward to derive this affine approximation:

```python
def sbox_affine(i:tuple):
    return (
        (0 + i[0]*1 + i[1]*1 + i[2]*0),
        (0 + i[0]*0 + i[1]*0 + i[2]*1),
        (1 + i[0]*0 + i[1]*2 + i[2]*2)
    )

SBOX_AFFINE = tuple(tri_to_int(sbox_affine(int_to_tyt(a)[0])) for a in range(27))
```

Plotting `SBOX_AFFINE`'s differntial table yields:

<p align="center">
  <img  src="rsrc/sbox_affine_diff.png" alt="SBOX_AFFINE differential table">
</p>

Also `SBOX` and `SBOX_AFFINE` turn out to be _really_ similay:

```python
print(SBOX_AFFINE)
print(SBOX)

# > (9, 10, 11, 1, 2, 0, 20, 18, 19, 3, 4, 5, 22, 23, 21, 14, 12, 13, 24, 25, 26, 16, 17, 15, 8, 6, 7)
# > (9, 10, 11, 1, 2, 0, 20, 18, 19, 3, 4, 5, 22, 23, 21, 14, 12, 26, 24, 25, 13, 16, 17, 15, 8, 6, 7)
```

`SBOX` and `SBOX_AFFINE` differ only in 2 values. This means that `25/27 ~ 93%` of the time, `SBOX` is affine!

Remember how the `SBOX` is the all important component that ensures that the cipher isn't affine? Having `SBOX` behave like an affine transform _most_ of the time is pretty bad.

Now, being affine `93%` of the time is very unlikely to make a particular instance of A3S encryption affine. Since the SBOX is applied _many many times_ during encryption, there is bound to be once where the SBOX isn't affine. 

So instead of saying that with this `SBOX`, the A3S approximates an affine transform, I believe it'll be way more helpful to know that small sections of an A3S encryption instance is gonna be affine. We don't know where and when these affine portions are, just that small parts of it is _pretty likely_ to be affine.

# Attack Approach

Unlike in the previous challenge [A3S](../A3S), it's unlikely that there's a plaintext-ciphertext pair that's completely affine for us to exploit. However, we can exploit the _very probably_ affine-ness of small parts of the encryption. 

Furthermore, we'd want to recover all the round keys, which are `3*8*9 = 216` `GF(3)` values. This means that we'd want to recover at least `216` constraints to solve for the round keys. Preferably, these constraints are linear to make solving for the round keys really easy.

There are 2 places to get these constraints. The obvious one is the key expansion algorithm, and the next is the plaintext-ciphertext pairs.

From now on, `kr` will refer to the round keys, plaintext is `pt` and ciphertext is `ct`.

So here's the plan:

1. Create a completely symbolic affine version of the A3S (that's a close approximation)
2. Recover as much constraints on `kr` from the key expansion.
3. Correlate the output of this affine version with the known `pt` and `ct` pairs to form more constraints on `kr`
4. With hopefully > 216 constraints, we can solve for `kr`.

Recovering all the round keys is equivalent to recovering the original key as the original key is simply the first few round keys.

## Creating a symbolic affine A3S

Just like in the [previous challenge](../A3S), I chose to modify the given implementation to accept `GF(3)` symbolic variables.

First, importing and processing what's given to us from the challenge:

```python
import mini as orig # Import challenge code
from mini import *  # Import challenge code
from out import *   # Import output of challenge code given to us

F = GF(3) # Field we are working in

# Convert known ct to the 3x3 tryte state
server_ct = [up(int_to_tyt(byt_to_int(ct)), W_SIZE ** 2, int_to_tyt(0)[0])[-1] for ct in all_enc]
# Flatten the ct to 27 trit vectors
server_ctv = [vector(F, [i for j in ct for i in j]) for ct in server_ct]
# Convert known pt to the 3x3 tryte state
server_pt = [up(int_to_tyt(pt), W_SIZE ** 2, int_to_tyt(0)[0])[-1] for pt in range(3^9)]
# Flatten the pt to 27 trit vectors
server_ptv = [vector(F, [i for j in pt for i in j]) for pt in server_pt]
```

Creating the symbolic variables:

```python
ptkr = PolynomialRing(F, ['pt%d'%i for i in range(9*3)] + ['kr%d'%i for i in range(9*3*8)]).gens()

# Flattened pt variables
pt_v = ptkr[:9*3]
# pt variables in 3x3 trit format
pt = tuple(tuple(pt_v[i*3:i*3+3]) for i in range(9))
# Flattened kr variables
kr_v = ptkr[9*3:]
# kr variables in 3x3 trit round keys format
kr = tuple(tuple(tuple(kr_v[i*9*3:i*9*3+9*3][j*3:j*3+3]) for j in range(9)) for i in range(8))
```

Now to modify the implementation to be `GF(3)` friendly:

```python
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

def unsbox_affine(i:tuple):
    return (
        (2 + i[0]*1 + i[1]*1 + i[2]*1),
        (1 + i[0]*0 + i[1]*2 + i[2]*2),
        (0 + i[0]*0 + i[1]*1 + i[2]*0)
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
```

Now I can simply call, say, `mix(pt)` and get the symbolic representation of `mix`:

<p align="center">
  <img  src="rsrc/symbolic_demo.png" alt="Symbolic model demo">
</p>

## Savaging ~168 linear constraints from the key expansion

Refering to the [.pdf](./chal/help.pdf) given:

<p align="center">
  <img  src="rsrc/key_expansion_pdf.png" alt="Key expansion algorithm">
</p>

The third line is an obvious linear constraint on `kr`. This generates `135` constraints:

```python
# xkey_mat1 * kr = xkey_const1

xkey_mat1 = []   # An array of 135 x 216 numbers
xkey_const1 = [] # A vector of 135 numbers

for i in range(6,24):
    if i%5 == 0: continue
    for j in range(9):
        r = vector([0]*3*8*9)
        r[i*9+j] = 1; r[i*9-9+j] = 2; r[i*9-9*5+j] = 2
        xkey_mat1.append(r)
        xkey_const1.append(0)
```

The second line is a little trickier. Notice the `sub` operation? It's affine `93%` of the time, but not _all_ the time. During the key expansion of the key for this challenge, `SBOX` is used a total of 12 times. This means there's around `12 * 8% ~ 1` instance where the constraint isn't linear.

If we replace the `sub` operation with its affine counterpart, we'd generate around `36` constraints. However, we expect about one of the SBOX applications to differ from our affine approximation. Each time that happens, it invalidates `3` of our constraints. So we should expect about `33` valid linear constraints, we just don't know _which_ of the SBOX affine approximation is invalid. No matter, since there's only `12` SBOX applications, we can simply bruteforce which one wasn't valid if it exists.

This gives a total of about `135 + 33 = 168` linear constraints on `kr`, still about `48` constraints short of the `216` goal.

## Scavaging the known pt-ct pairs for 54 linear constraints. 

This is where the bulk of my time was spent solving this challenge: Scraping for enough constraints to proceed.

We'll start by defining the areas of interest. These areas are where we'd correlate the output of our affine model and the actual A3S output to derive more constraints. I've found 2 areas of interests, each giving `27` constraints. There are more areas, but the areas I've chosen make implementation really easy:

```python
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
```

As you can see, I chose to stop at checkpoints `0` and `1`, where `0` is right before the last `sub` operations, and `1` is before the last `add_round_key` operation.

I've also implemented a function that outputs matrices of interest at that checkpoint:

```python
def gen_mat(checkpoint):

    """
    Generates matrices and vectors corresponding
    to the affine transform to `checkpoint`.

    Returns ptmat, ptconst, kptmat where
        ct = ptmat*pt + ptconst + kptmat*kr
    """
    
    ctt = forward_to_checkpoint(pt, kr, checkpoint)
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
```

### Checkpoint 0

