# The Challenge

You are given `code.py` which encrypts some text with ECC, `enc.txt` which presumably contains the encrypted text. In addition, there is `partial_plaintext.txt` which contains known plaintext from index 210 to 242 and an image `partial_plaintext.png` that contains more segments of known plaintext. I ended up not using `partial_plaintext.png` so I'm not sure why it was included. These files are included in the folder [Challenge Files](https://github.com/JuliaPoo/Collection-of-CTF-Writeups/edit/master/CDDC2020/A%20Kind%20of%20Crypto/Challenge%20Files).

Presumably the goal is to decrypt `enc.txt`, which is encrypted with a different seed than given in `code.py`.

## Solution

There are two points `Q` and `P` defined on a Elliptic Curve `E`. From a seed `s`, it computes `r = sQ` and a new seed `sP` to be used in the next iteration. The sequence of `r`s generated are used to encrypt the plaintext by bit-wise xoring. However, each time only the `30 * 8` LSB of `r` are used to decrypt a 30 byte chunk of plaintext each time. `partial_plaintext.txt` hence contain plaintext of the 7th chunk and the first 3 bytes of the 8th chunk.

We can recover `r` from `partial_plaintext.txt` with ease simply by bit-wise xoring. However, there are 16 bits of `r` that seem to be lost during the encryption process. Fortunately `2**16 = 65536` is a pretty small number to bruteforce so we can recover the whole `r` no problem.

Recovering `s` via `r = sQ` is a little more challenging. In a proper implementation this would be impossible except via bruteforce due to the [Discrete Logarithm Problem](https://en.wikipedia.org/wiki/Elliptic-curve_cryptography). This turns out to be impractical with this challenge (I tried). However, in the event that the order of the curve is equal to the order of the finite field, i.e. equal to the prime used in `code.py`, there is something called the [Smart's Attack](https://wstein.org/edu/2010/414/projects/novotney.pdf) that enables recovering of `s` from `r = sQ` in linear time. An implementation of Smart's Attack can be found in the Python library [ecpy](https://github.com/elliptic-shiho/ecpy).

Now, we can directly apply Smart's Attack to recover `s`. However, remember that we have to bruteforce a maximum of 65536 different `r`. It proves to be slow to apply Smart's Attack that many time. We can get around this by computing a value `e` where `eQ = P`, which requires one application of Smart's Attack, and then `s'` for the next iteration can simply be calculated as `s' = sP = s * eQ = e * sQ = e * r`. We can then use this new `s'` to decrypt the first 3 bytes of the 8th chunk and compare it with that in `partial_plaintext.txt`.

```python
from ecpy import * # https://github.com/elliptic-shiho/ecpy

p = 112817876910624391112586233842848268584935393852332056135638763933471640076719

A = 49606376303929463253586154769489869489108883753251757521607397128446713725753
B = 79746959374671415610195463996521688925529471350164217787900499181173830926217
P = 112817876910624391112586233842848268584935393852332056135638763933471640076719
Px = 103039657693294116462834651854367833897272806854412839639851017006923575559024
Py = 77619251402197618012332577948300478225863306465872072566919796455982120391100
Qx = 54754931428196528902595765731417656438047316294230479980073352787194748472682
Qy = 31061354882773147087028928252065932953521048346447896605357202055562579555845
F = FiniteField(p)
E = EllipticCurve(F, A, B)

Q = E(Qx, Qy)
P = E(Px, Py)

e = SSSA_Attack(F, E, Q, P)
assert Q * e == P
print ("[+] e = %d" % e)

enc = open('Challenge Files/enc.txt', 'r').read()
ct = bytearray([int(enc[2*i:2*i+2], 16) for i in range(len(enc)//2)])[210:] # Get bytearray from enc with offset 210
pt = b'been done and we are seeing extre'

pt30 = pt[:30]
ct30 = ct[:30]
r30 = [x^y for x,y in zip(pt30, ct30)]

def get_r(pt30, ct30, r30):

    'Generates all possible values of r'
    
    m = int(r30[0])
    for c in r30[1:]:
        m *= 256
        m += c
      
    x = m
    i = 0
    while True:
        i += 1
        x += 2**(30*8) #* 60359
        y = E.get_corresponding_y(x)
        if y:
            yield i, E(x,y)

for i, r in get_r(pt30, ct30, r30):
    print(i, end='\r')
    s2 = (e * r).x.x
    r2 = (Q * s2).x.x & (2**(30*8) - 1)
    pt_dec = []
    for j in range(1,4):
        pt_dec += bytearray([((r2>>((30-j)*8))&0xff)^ct[30*1+j-1]])
    pt_dec = bytes(pt_dec)
        
    if pt_dec == pt[30:]:
        break

# >>> [+] e = 72529124805871229360171330254593943220566431521453438361067644203504289580075
# >>> 60359
```

Alright! So we've found both `r` and `s'`! We can simply feed `s'` as a seed into the decryption function in `code.py` to decrypt everything from the 9th chunk onwards!

```python
from code import *

ec = EC(A,B,p)
P0 = Point(Px, Py)
Q0 = Point(Qx, Qy)

s = STREAM(ec, s2, P0, Q0)
print(s.decryption(ct[60:]).decode('utf-8'))

# up to 98%. We are ready to inject this odourless, colourless and tasteless liquid into all our water pumps. Prepare yourselves for CDDC20{maS5_brA1nwashINg_anD_w0rLD_dOMINA7ioN}!! HAHAHAHAHHAAAAA cheers to the success of our evil planz!!!
```

## Beyond: Recovering the whole message and the original seed used

I wondered whether there could be some funny secret message the author placed in the other parts of the encrypted message, and maybe in the seed? Good thing is, from here, all we gotta do is apply Smart's Attack 9 more times to recover the original seed. However, at every iteration, there are two solutions for `s`, namely `s` and `p-s`, where `p` is the modulo of the finite field, which isn't a problem, just bruteforce all solutions.

```python
s1 = s2
s1 =     SSSA_Attack(F, E, P, E(s1, E.get_corresponding_y(s1)))
s1 =     SSSA_Attack(F, E, P, E(s1, E.get_corresponding_y(s1)))
s1 = p - SSSA_Attack(F, E, P, E(s1, E.get_corresponding_y(s1)))
s1 =     SSSA_Attack(F, E, P, E(s1, E.get_corresponding_y(s1)))
s1 = p - SSSA_Attack(F, E, P, E(s1, E.get_corresponding_y(s1)))
s1 = p - SSSA_Attack(F, E, P, E(s1, E.get_corresponding_y(s1)))
s1 = p - SSSA_Attack(F, E, P, E(s1, E.get_corresponding_y(s1)))
s1 = p - SSSA_Attack(F, E, P, E(s1, E.get_corresponding_y(s1)))
s1 = p - SSSA_Attack(F, E, P, E(s1, E.get_corresponding_y(s1)))

ct = bytearray([int(enc[2*i:2*i+2], 16) for i in range(len(enc)//2)])

s = STREAM(ec, s1, P0, Q0)
print(s.decryption(ct).decode('utf-8'))

# >>> Greetingzz, my fellow comrades from Unduplicitous Corp. Our researchers have made tremendous progress with Brainwashing Agent 2.0 - it is now alot more effective than the previous version. Thorough testing has been done and we are seeing extremely high success rates of up to 98%. We are ready to inject this odourless, colourless and tasteless liquid into all our water pumps. Prepare yourselves for CDDC20{maS5_brA1nwashINg_anD_w0rLD_dOMINA7ioN}!! HAHAHAHAHHAAAAA cheers to the success of our evil planz!!!

print('Possible seed used:')
print('\tdec: %d'%(s1))
print('\thex: %s'%hex(s1)[2:])

# >>> Possible seed used:
# >>>	dec: 151386116941948313024325411690796683746
# >>>	hex: 71e3e7d3fac11617c282d57c4ab211e2
```

So unfortunately no secret message in the seed but it was fun.
