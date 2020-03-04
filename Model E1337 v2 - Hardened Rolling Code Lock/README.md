# Model E1337 v2 - Hardened Rolling Code Lock - FLAG0

The solution for this flag is very similar to [Model E1337 - Rolling Code Lock](https://github.com/JuliaPoo/Hacker101-CTF/blob/master/Model%20E1337%20-%20Rolling%20Code%20Lock/README.md). Hence this writeup will be less detailed.

## Preamble

the source code of the PRNG can be found by appending ```/rng``` to the URL.
Just like Model E1337 - Rolling Code Lock, the ```getrandbit``` function can be defined as

```python
def getrandbit():
    global state
    
    ret <<= 1
    ret |= state & 1
    
    for k in range(3):
        state = (state << 1) ^ (state >> 61)
        state &= 0xFFFFFFFFFFFFFFFF
        state ^= 0xFFFFFFFFFFFFFFFF

        for j in xrange(0, 64, 4):
            cur = (state >> j) & 0xF
            cur = (cur >> 3) | ((cur >> 2) & 2) | ((cur << 3) & 8) | ((cur << 2) & 4)
            state ^= cur << j
```

## Reversing the PRNG

Tracking the bits of the state is done similarly to Model E1337 - Rolling Code Lock

```python

# This bit index represents '1'
RESERVED = {64}

def initialise_state_a():
    state_a = [{63-i} for i in range(64)]
    return state_a

def getrandbit_track2(state_a):
    
    '''
    Same as getrandbit_track but for the 2nd challenge
    '''
    
    ret = state_a[-1]

    for k in range(3):
        new = [None] * 64
        for i in range(61):
            new[i] = state_a[i + 1]
        new[61] = state_a[62] ^ state_a[0]
        new[62] = state_a[63] ^ state_a[1]
        new[63] = state_a[2]
        state_a = new[:]

        state_a = [i ^ RESERVED for i in state_a]

        for j in range(0, 64, 4):
            if j == 0: cur = state_a[-j-4:]
            else: cur = state_a[-j-4:-j]
            cur_ = [cur[-1], cur[-1], cur[-4], cur[-4]]
            i = 0
            for k in range(60-j, 64-j):
                state_a[k] = state_a[k] ^ cur_[i]
                i += 1
    
    return state_a, ret
    
```

## Predicting

The rest is the same as my solution for Model E1337 - Rolling Code Lock

```python
import time
t = time.time()

rng_res = [2288066835647663286]  # numbers from the server
NBIT = 64                        # Number of bits of each number
NNUMBER = len(rng_res)           # Number of numbers given
    
# Create equations that corresponds to rng_res
ret_b = [int(i) for i in "".join(['0'*(NBIT-len(bin(i)[2:])) + bin(i)[2:] for i in rng_res])]
all_ret = []
state_a = initialise_state_a()
for i in range(NBIT*NNUMBER):
    state_a, ret = getrandbit_track2(state_a)
    all_ret.append(ret)
    
eqns = [[x,y] for x,y in zip(all_ret, ret_b)]
eqns = solve_xor(eqns)
eqns_rhs = [eq[1] for eq in eqns]

# Create expressions for the bits of the next number
pred_eqns = []
for i in range(NBIT):
    state_a, ret = getrandbit_track2(state_a)
    pred_eqns.append(ret)

# decompose pred_eqns to the known eqns (from server)
# to predict next number
bits = ""
for pred_eqn in pred_eqns:
    composition = decompose_eqn(eqns, pred_eqn)
    pred = eqns_rhs[composition[0]]
    for idx in composition[1:]:
        pred = pred ^ eqns_rhs[idx]
    bits = bits + str(pred)
    
print("Predicted:", int(bits, 2))
print("Time taken:", "{}s".format(round(time.time() - t,2)))
```

Output:
```
Predicted: 10951372415822408502
Time taken: 0.15s
```

The full code is available in ```Solve.py```
