import dlp
import time
from sympy.ntheory import randprime
import gmpy2
from gmpy2 import mpz
# mpz_urandomb, random_state

random_state = gmpy2.random_state(time.time_ns())
# p = int(input("Введите p: "))
p = randprime(2**40, 2**60)
# g = int(input("Введите g: "))
g = gmpy2.mpz_urandomb(random_state, 60)
# h = int(input("Введите h: "))
h = gmpy2.mpz_urandomb(random_state, 60)
print("Solving g ^ x = h (mod p), where x = log_g h (mod N), ord(g) = N")
print(g, "^ x =", h, "( mod", p, ")")
print("p length (bits):", p.bit_length())
print("g length (bits):", g.bit_length())
print("h length (bits):", h.bit_length())
ord = dlp.order(g, p)
start = time.perf_counter()
x = dlp.pohlig_hellman(p, g, h, ord)
stop = time.perf_counter()
print("\nExecution time ", stop - start, "s")
if x == -1 and ord != p - 1:
    print("No solution")
else:
    print("x =", x, "mod", ord)
    if pow(g, x, p) != h % p: # and ord == p - 1:
        print("Failed")
    else:
       print("Passed")
