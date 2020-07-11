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
    
import time
t = time.time()

rng_res = [63235802, 41186890]   # numbers from the server
NBIT = 26                        # Number of bits of each number
NNUMBER = len(rng_res)           # Number of numbers given
    
# Create equations that corresponds to rng_res
ret_b = [int(i) for i in "".join(['0'*(NBIT-len(bin(i)[2:])) + bin(i)[2:] for i in rng_res])]
all_ret = []
state_a = initialise_state_a()
for i in range(NBIT*NNUMBER):
    state_a, ret = getrandbit_track(state_a)
    all_ret.append(ret)
    
eqns = [[x,y] for x,y in zip(all_ret, ret_b)]
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