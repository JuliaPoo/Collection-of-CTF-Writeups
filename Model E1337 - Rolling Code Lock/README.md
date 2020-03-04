# Model E1337 - Rolling Code Lock - FLAG0

Easy to find the solution to this flag online, so I'm not gonna repeat the work.

# Model E1337 - Rolling Code Lock - FLAG1

## Preamble

The source code for the PRNG is

```python
def setup(seed):
    global state
    state = 0
    for i in xrange(16):
        cur = seed & 3
        seed >>= 2
        state = (state << 4) | ((state & 3) ^ cur)
        state |= cur << 2

def next(bits):
    global state

    ret = 0
    for i in xrange(bits):
        ret <<= 1
        ret |= state & 1
        state = (state << 1) ^ (state >> 61)
        state &= 0xFFFFFFFFFFFFFFFF
        state ^= 0xFFFFFFFFFFFFFFFF

        for j in xrange(0, 64, 4):
            cur = (state >> j) & 0xF
            cur = (cur >> 3) | ((cur >> 2) & 2) | ((cur << 3) & 8) | ((cur << 2) & 4)
            state ^= cur << j

    return ret
```

The server calls ```next(26)``` to generate the next key. Defining a function ```getrandbit()```

```python
def getrandbit():
    global state
    
    ret <<= 1
    ret |= state & 1
    state = (state << 1) ^ (state >> 61)
    state &= 0xFFFFFFFFFFFFFFFF
    state ^= 0xFFFFFFFFFFFFFFFF

    for j in xrange(0, 64, 4):
        cur = (state >> j) & 0xF
        cur = (cur >> 3) | ((cur >> 2) & 2) | ((cur << 3) & 8) | ((cur << 2) & 4)
        state ^= cur << j
```

Predicting the output of ```next(26)``` is equivalent to predicting the output of ```getrandbit()``` 26 times in a row.

## Reversing the PRNG

My approach to predict the ```getrandbit``` is to use the output of a few ```next(26)``` to recover as much of ```state``` as needed to be able to predict the next bits. An observation is that the output of ```getrandbit``` is can be represented as a series of XOR of the bits of the original state (from ```setup(seed)```).

Define the original state as a list of sets. ```state_a = [{63}, {62}, {61}, {60}...{0}]```. The numbers represent the index of the bit of the state, ```63``` being the ```MSB``` and ```0``` being the ```LSB```. As ```state_a``` gets put through ```getrandbit```, it will evolve. E.g. ```state_a = [{0}, {63, 61}...]```. The nth set in ```state_a``` represents the nth bit of ```state_a``` as a sequence of XORs of the bits of the original state. E.g. if ```state_a[-5] = {63, 61}```, then ```state & 2**4 = a63 ^ a61```, where the original state has bits ```a63, a62, a61..., a0```.

A function ```getrandbit_track(state_a)``` would evolve ```state_a``` upon passing through ```getrandbit```. My implementation as follows. Do note that I reserved index ```64``` to represent ```1```, while all the other index represents the nth bit of the original state:

```python
# This bit index represents '1'
RESERVED = {64}

def initialise_state_a():
    state_a = [{63-i} for i in range(64)]
    return state_a
    
def getrandbit_track(state_a):
    
    '''
    Input:
    state_a
        A representation of the prng state
        state_a is a list of sets:
            The nth set represents the nth bit of the state
            Each set is made of index of the original state
            
        E.g. state_a[3] = {1,4,5,63}:
        If the original state has bits (MSB) a63, a62, a61, a60..., a0 (LSB)
        state & 2**2 = a1 ^ a4 ^ a5 ^ a63
        
    Output:
    returns state_a that represents the state after generating a random bit
    '''
    
    # ret = state & 1
    ret = state_a[-1]

    # state = (state << 1) ^ (state >> 61)
    ##########################################################################
    # a63, a62 ... a4,  a3,  a2,  a1,  a0  # state
    # a62, a61 ... a3,  a2,  a1,  a0,   0  # state <, 1
    #   0,   0 ...  0,   0, a63, a62, a61  # state >> 61
    # a62, a61 ... a3,  a2,  b1,  b2, a61  # (state << 1) ^ (state >> 61) = C
    # a1 = b1 ^ a63
    # a0 = b2 ^ a62
    ##########################################################################
    new = [None] * 64
    for i in range(61):
        new[i] = state_a[i + 1]
    new[61] = state_a[62] ^ state_a[0]
    new[62] = state_a[63] ^ state_a[1]
    new[63] = state_a[2]
    state_a = new[:]

    # state ^= 0xFFFFFFFFFFFFFFFF
    state_a = [i ^ RESERVED for i in state_a]

    # for j in range(0, 64, 4):
    #         cur = (state >> j) & 0xF
    #         cur = (cur >> 3) | ((cur >> 2) & 2) | ((cur << 3) & 8) | ((cur << 2) & 4)
    #         state ^= cur << j
    ##########################################################################
    # 0 ... c3 c2 c1 c0             # cur
    # 0 0 0  0 ...  0  0  0  0 c3   # (cur >> 3)
    # 0 0 0  0 ...  0  0  0 c3  0   # ((cur >> 2) & 2)
    # 0 0 0  0 ...  0 c0  0  0  0   # ((cur << 3) & 8)
    # 0 0 0  0 ...  0  0 c0  0  0   # ((cur << 2) & 4)
    # 0 0 0  0 ...  0 c0 c0 c3 c3   # (cur >> 3) | ((cur >> 2) & 2) | ((cur << 3) & 8) | ((cur << 2) & 4)
    ##########################################################################
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

Using the evolved  ```state_a```, we can represent the result of ```getrandbit``` in terms of the bits of the original state. By getting numbers from the server, we can have a system of equations to solve for the original state. Since the only operation in this system of equations is XOR, we can use Gaussian Elimination to solve. My implementation as follows:

```python
def solve_xor(eqns):
    
    '''
    Input: 
    eqns is a list of eqn:
        Each eqn is [set of bit index, 1 or 0]
        E.g. [{0,1,4,5,2}, 1]
        Which corresponds to:
            a0 ^ a1 ^ a4 ^ a5 ^ a2 = 1
            
    Output:
    Goes Gaussian Elimination of the eqns and output final system of eqns
    Should the system of eqns not be valid, a warning is given but execution does not stop
    '''
    
    # Remove RESERVED
    res = list(RESERVED)[0]
    for idx, eqn in enumerate(eqns):
        if res in eqn[0]: 
            eqn[0] = eqn[0] ^ RESERVED
            eqn[1] = eqn[1] ^ 1
            eqns[idx] = eqn

    # Solve:
    used = []
    for bit in range(64):

        # Find first equation with bit
        cut_short = False
        eqn_idx = 0

        while eqn_idx < len(eqns):
            eqn = eqns[eqn_idx]
            if bit not in eqn[0]:
                eqn_idx += 1
                continue
            elif eqn_idx in used:
                eqn_idx += 1
                continue
            else: 
                first_idx = eqn_idx
                eqn_idx += 1
                break
            if eqn_idx == len(eqns) - 1:
                cut_short = True

        if cut_short: continue

        used.append(first_idx)

        # Eliminate bit in rest of equations
        eqn_idx = 0
        while eqn_idx < len(eqns):
            eqn = eqns[eqn_idx]
            if first_idx == eqn_idx: 
                eqn_idx += 1
                continue
            if bit in eqn[0]:
                mod_eqn = eqns[first_idx][0] ^ eqn[0]
                mod_bit = eqns[first_idx][1] ^ eqn[1]
                eqns[eqn_idx] = [mod_eqn, mod_bit]
            eqn_idx += 1

    # Check if system of equations is consistent
    for eqn in eqns:
        if len(eqn[0]) == 0 and eqn[1] == 1:
            print("[WARNING] System of equations make no sense")

    # Remove empty equations
    eqns = [i for i in eqns if len(i[0]) != 0]
    
    return eqns
 ```
 
 However, it seems to be impossible to recover all the bits of the original state. It seems that after a certain point, any additional number given by the server will not solve more bits of the original state. An example of the final equations are represented below:
 
 ```python
a0  =  1
a2 ^ a61  =  1
a1 ^ a62 ^ a63  =  1
a59 ^ a61  =  1
a55 ^ a57  =  1
a51 ^ a53  =  1
a47 ^ a49  =  0
a43 ^ a45  =  0
a41 ^ a39  =  1
a35 ^ a37  =  0
a33 ^ a31  =  0
a27 ^ a29  =  1
a25 ^ a23  =  0
a19 ^ a21  =  0
a17 ^ a15  =  0
a11 ^ a13  =  0
a9 ^ a7  =  0
a3 ^ a5  =  1
a5 ^ a6  =  0
a4 ^ a5  =  0
a8 ^ a9  =  0
a9 ^ a10  =  0
a13 ^ a12  =  0
a13 ^ a14  =  1
a17 ^ a18  =  0
a16 ^ a17  =  1
a20 ^ a21  =  1
a21 ^ a22  =  1
a24 ^ a25  =  1
a25 ^ a26  =  1
a29 ^ a30  =  1
a28 ^ a29  =  0
a32 ^ a33  =  1
a33 ^ a34  =  1
a36 ^ a37  =  1
a37 ^ a38  =  1
a41 ^ a42  =  0
a40 ^ a41  =  0
a44 ^ a45  =  0
a45 ^ a46  =  1
a48 ^ a49  =  0
a49 ^ a50  =  0
a53 ^ a54  =  0
a52 ^ a53  =  1
a56 ^ a57  =  0
a57 ^ a58  =  0
a60 ^ a61  =  0
a61 ^ a62  =  0
```

## Predicting

However, it is not needed to recover all the bits of the original state to predict ```getrandbit```. Since an additional number from the server does not give any more information about the original bits (in other words, does not add anything to the above system of equations), it means that the value of the next number can be completely determined by the above system of equations. In particular, given an expression that represents the result of ```getrandbit```, it should be possible to represent said expression as a series of XORs from the equations given above, and therefore be able to predict the result of ```getrandbit```. I achieved this with gaussian elimination again, shown below:

```python
def decompose_eqn(eqns_, eqn_to_decompose):
    
    '''
    Input:
    eqns_: a list of eqn
    eqn_to_decompose: Same as an eqn but without the result. E.g. {1,5,3,2}
    
    Output:
    Tries to generate a list of index of eqns from eqns_ which when xored 
    together produces eqn_to_decompose
    '''
    
    A = eqn_to_decompose
    contains_res = RESERVED in A
    
    eqns = eqns_ + [[A]]
    
    used = []
    bits_A = list(A)

    composition = []
    for bit in range(64):
        
        # Find first equation with bit
        cut_short = False
        eqn_idx = 0
        
        while eqn_idx < len(eqns)-1:
            eqn = eqns[eqn_idx]
            if bit not in eqn[0]:
                eqn_idx += 1
                continue
            elif eqn_idx in used:
                eqn_idx += 1
                continue
            else: 
                first_idx = eqn_idx
                break
            if eqn_idx == len(eqns) - 2:
                cut_short = True
                
        if cut_short: continue

        used.append(first_idx)
        
        # Gaussian eliminate only the last row (eqn_to_decompose)
        if bit in eqns[-1][0]:
            mod_eqn = eqns[first_idx][0] ^ eqns[-1][0]
            eqns[-1] = [mod_eqn, None]
            composition.append(first_idx)

    # Check if equation can be composed off the other eqns
    if len(eqns[-1][0]) != 0: 
        print("[WARNING] Equation cannot be composed from given eqns")
        
    return composition
 ```
 
 Putting it all together, here is my code for predicting ```next(26)``` given ```2``` numbers from the server:
 
 ```python
 import time
t = time.time()

rng_res = [63235802, 41186890]   # numbers from the server
NBIT = 26                        # Number of bits of each number
NNUMBER = len(rng_res)           # Number of numbers given
    
# Create equations that corresponds to rng_res
ret_b = [int(i) for i in "".join(['0'*(NBIT-len(bin(i)[2:])) + bin(i)[2:] for i in rng_res])]
ALL_RET = []
state_a = initialise_state_a()
for i in range(NBIT*NNUMBER):
    state_a, ret = getrandbit_track(state_a)
    ALL_RET.append(ret)
    
eqns = [[x,y] for x,y in zip(ALL_RET, ret_b)]
eqns = solve_xor(eqns)
eqns_rhs = [eq[1] for eq in eqns]

# Create expressions for the bits of the next number
pred_eqns = []
for i in range(NBIT):
    state_a, ret = getrandbit_track(state_a)
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
Predicted: 1352947
Time taken: 0.05s
```

The whole python code can be found in ```Solve.py```
