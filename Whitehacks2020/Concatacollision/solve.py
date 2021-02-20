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