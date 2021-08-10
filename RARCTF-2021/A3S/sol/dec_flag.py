from a3s import *

# Key recovered from key_recover.sage
key_rec = ((2, 1, 1), (1, 1, 1), (1, 1, 2), (1, 0, 0), (2, 2, 0), (2, 2, 1), (0, 1, 2), (2, 1, 2), (0, 2, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0))

# Modify d_a3s to use key in tyt form
def d_a3s(ctt, key):
    c       = byt_to_int(ctt)
    c       = up(int_to_tyt(c), W_SIZE ** 2, int_to_tyt(0)[0])[-1] # Fixed block size
    keys    = expand(key)[::-1] # tryte array
    assert len(keys) == KEYLEN

    msg = c
    msg = T_uxor(msg, keys[0])

    for r in range(1, len(keys) - 1):
        msg = tuple([msg[i] for i in UN_SHIFT_ROWS])    # UN SHIFT...
        msg = u_sub_wrd(msg)                            # UN SUB...
        msg = T_uxor(msg, keys[r])                      # UN ADD...
        msg = mix_columns(msg, I_CONS)                  # UN MIX!

    msg  = tuple([msg[i] for i in UN_SHIFT_ROWS])
    msg  = u_sub_wrd(msg)
    msg  = T_uxor(msg, keys[-1])                     # last key

    msg = tyt_to_int(msg)
    return int_to_byt(msg)

flag_enc = b'\x01\x00\xc9\xe9m=\r\x07x\x04\xab\xd3]\xd3\xcd\x1a\x8e\xaa\x87;<\xf1[\xb8\xe0%\xec\xdb*D\xeb\x10\t\xa0\xb9.\x1az\xf0%\xdc\x16z\x12$0\x17\x8d1'

out = []
for i in chunk(flag_enc):
    out.append(d_a3s(i, key_rec))
    
print("Flag:", unchunk(out).decode())