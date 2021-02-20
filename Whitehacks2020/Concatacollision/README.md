# The Challenge

You are given a file [concat.py](./dist/concat.py) that's running on their server. It randomly generates two numbers (`x0` and `x1`) and asks you to input a third number (`s`) which is validated via the `check` function.

## Metadata

I did not participate in Whitehacks 2020. I was bored and did a challenge. The reason why I'm writing this is because the [canon solution given by the challenge authors](https://github.com/Whitehat-Society/whitehacks-challenges-2020-public/blob/master/challenges/crypto/Concatacollision/token.py) is super clunky. 

## Solution

The function `check(s:int, x0:int, x1;int)` returns true if `int(str(int(str(x1) + str(s)) % x0) + str(s)) % x1 == 1` (I modified the argument types so all inputs are integers). This can be formalised as:

```
x1 * 10^p + s == a  mod x0
a  * 10^p + s == 1  mod x1
```

Where `p == len(str(s))`, or the number of digits `s` has. Let's combine this system into a single equation via Chinese Remainder Theorem:

```
s == a - x1 * 10^p  mod x0
s == 1 - a  * 10^p  mod x1

⇓

s == (1 - a*10^p)*modinv(x0,x1)*x0 + (a - x1*10^p)*modinv(x1,x0)*x1
  == a * (modinv(x1,x0)*x1 - modinv(x0,x1)*x0*10^p) 
    + (modinv(x0,x1)*x0 - modinv(x1,x0)*x1*x1*10^p)
  == a * N + K  mod x0*x1
```

Where `modinv(a,b)` is the modular inverse of `a` in `mod b`. Now, this obviously only works when the modular inverses exists, which requires `gcd(x0,x1)==1`. However, we only need it to work _sometimes_. It if doesn't work, you can simple restart the connection and get a new code.

As above, every value of `a` will give a solution for `s`. So why not set `a==0` to simplify things?

```
s == K
  == modinv(x0,x1)*x0 - modinv(x1,x0)*x1*x1*10^p mod x0*x1
```

However, there is still one unknown variable here that's dependent on `s`, which is `p`, the number of digits of `s`. Since we only need it to work _sometimes_, why not set `p==len(str(x0*x1))` as that happens fairly often.

Hence the (super duper short) code to solve for `s`:

```py
# Python 3.8 or above

x0 = 252407470330613386741277095757154302720 # from server
x1 = 177907764183659887207481969886517086087 # from server
check = lambda s: int(str(int(str(x1) + str(s)) % x0) + str(s)) % x1

p = len(str(x0*x1))
s = pow(x0,-1,x1)*x0 - pow(x1,-1,x0)*x1*x1*10**p
s = s % (x0*x1)

assert check(s)==1, len(s) == p-1
print(s)

# > 31175392867133471711106087115764756200474617904802451549159453225242446295040
```

Input this number into the server and we get the flag.
