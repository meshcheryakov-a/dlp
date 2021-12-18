import bts
from sympy.ntheory import primefactors, factorint
import gmpy2
from gmpy2 import mpz
# add, divexact, f_mod, mpz_urandomb + random_state, mul, powmod, isqrt

# вычисляет порядок элемента g в группе F*_p
def order(g, p):
    # полагаем, что порядок максимальный
    t = p - 1
    while 1:
        # разложение на простые множители числа t (каноническое разложение)
        factors = primefactors(t)
        ord = t
        for i in range(len(factors)):
            # t_i делит p-1
            # t_i = t // factors[i]
            t_i = gmpy2.divexact(mpz(t), mpz(factors[i]))
            # g_i = pow(g, t_i, p)
            g_i = gmpy2.powmod(mpz(g), mpz(t_i), mpz(p))
            # если g^t_i сравнимо с 1, то порядок g можно уменьшить до t_i
            if g_i == 1:
                t = t_i
                break
        # в случае, если g^t_i не сравнимо с 1 для всех t_i делителей p-1, то порядок g равен t
        if ord == t:
            return ord

# решение сравнения g^x=h(mod p) методом BSGS, ord(g)=N
def shanks(p, g, h, N):
    # выбор n
    # n = 1 + int(pow(N, 1/2))
    n = 1 + gmpy2.isqrt(mpz(N))
    # l = []
    # x=a+b*n
    a = 0
    b = 0
    # формирование списка (бинарного дерева поиска) BS, g^i, 0<=i<n
    g_i = 1
    l = bts.BTS(0, g_i)
    for i in range(n):
        # l.append(pow(g, i, p))
        # l.append(g_i)
        l.insert(i, g_i)
        # g_i = (g_i * g) % p
        g_i = gmpy2.f_mod(mpz(gmpy2.mul(mpz(g_i), mpz(g))), mpz(p))
    # поиск мультипликативного обратного для g и возведение его в степень n
    # g_inv = pow(g, -n, p)
    g_inv = gmpy2.powmod(mpz(g), mpz(-n), mpz(p))
    # формирование списка GS, h*g^-nj, 0<=j<n
    # g_0 = h % p
    # g_j = gmpy2.f_mod(mpz(h), mpz(p))
    g_j = h
    for j in range(n):
        # g_j = (h * pow(g_inv, j, p)) % p
        # проверка на пересечение со списком BS
        # if g_j in l:
        ind = l.search(g_j)
        if ind != -1:
            b = j
            # a = l.index(g_j)
            a = ind
            break
        # пересечение списков не найдено, логарифм не существует
        if ind == -1 and j == n - 1:
            return -1
        # g_j = (g_j * g_inv) % p
        g_j = gmpy2.f_mod(mpz(gmpy2.mul(mpz(g_j), mpz(g_inv))), mpz(p))
    return a + b * n

# алгоритм сведения к собственным подгруппам, решающий сравнение g^x=h(mod p), ord(g)=N=p'^e с помощью алгоритма Шенкса
def group(p, g, h, N):
    x = []
    # определяем p' и e, в данном случае ord=p'
    factors = factorint(N)
    ord, e = sorted(factors.items())[0]
    if e == 1:
        return shanks(p, g, h, N)
    # представляем x=x_0+x_1*p'+x_2*p'^2+...+x_(e-1)*p'^(e-1)
    # g_i = pow(g, N // ord, p)
    g_i = gmpy2.powmod(mpz(g), mpz(gmpy2.divexact(mpz(N), mpz(ord))), mpz(p))
    # проходим циклом по всем x_k, 0<=k<e, решая сравнение g_(k+1)^x_k=h_(k+1)(mod p) с помощью алгоритма Шенкса
    # g_(k+1)=g_1=g^(p'^(e-1)) для любых k, поэтому вычисляется до начала цикла
    for i in range(1, e + 1):
        # вычислям значение показателя p'^(e-i)
        exp = pow(ord, e - i)
        # h_i = pow(h, exp, p)
        h_i = gmpy2.powmod(mpz(h), mpz(exp), mpz(p))
        # формируем значение h_i, которое будет зависеть от произведений g^(-x_j*p'^(e-i+j)), 0<=j<k
        # например, h_1=h^(p^(e-1)), h_2=h^(p^(e-2))*g^(-x_0*p^(e-2)), h_3=h^(p^(e-3))*g^(-x_0*p^(e-3))*g^(-x_1*p^(e-2)) и т.д.
        for j in range(len(x)):
            exp = pow(ord, e - i + j)
            # h_i *= pow(g, -x[j] * exp, p)
            h_i = gmpy2.mul(mpz(h_i), mpz(gmpy2.powmod(mpz(g), mpz(-x[j] * exp), mpz(p))))
        # h_i = h_i % p
        h_i = gmpy2.f_mod(mpz(h_i), mpz(p))
        # решаем простое сравнение g_(k+1)^x_k=h_(k+1)(mod p) с помощью алгоритма Шенкса, ord(g_(k+1))=p' будет простым числом
        # откуда находим x_k
        x_k = shanks(p, g_i, h_i, ord)
        # если решения нет, то возвращаем -1
        if x_k == -1:
            return -1
        x.append(x_k)
    res = 0
    # восстанавливаем x=x_0+x_1*p'+x_2*p'^2+...+x_(e-1)*p'^(e-1)
    for j in range(len(x)):
        res += x[j] * pow(ord, j)
    return res

# решение сравнения g^x=h(mod p), где ord(g)=N=p_1^e_1 * p_2^e_2 * ... * p_m^e_m
# сведением исходной задачи к набору задач DLP с ord(g_i)=p^e
def pohlig_hellman(p, g, h, N):
    # факторизуем порядок ord(g)=N
    factors = factorint(N)
    print("ord(g) = N =", N, "=", factors)
    y = []
    # редукция
    # проходим по всем сомножителям p_i^e_i числа N, находим N_i=N/p_i^e_i, g_i=g^N_i, h_i=h^N_i
    # решаем сравнение g_i^y_i=h_i(mod p_i^e_i) с помощью алгоритма сведения к собственным подгруппам
    i = 1
    for p_i, e_i in sorted(factors.items()):
        i = str(i)
        print("\ni =", i)
        ord = pow(p_i, e_i)
        print("p_" + i, "=", p_i)
        print("e_" + i, "=", e_i)
        # N_i = N // ord
        N_i = gmpy2.divexact(mpz(N), mpz(ord))
        print("N_" + i, "=", N_i)
        # g_i = pow(g, N_i, p)
        g_i = gmpy2.powmod(mpz(g), mpz(N_i), mpz(p))
        print("g_" + i, "=", g_i)
        # h_i = pow(h, N_i, p)
        h_i = gmpy2.powmod(mpz(h), mpz(N_i), mpz(p))
        print("h_" + i, "=", h_i)
        print("Solving g_" + i, "^ y_" + i, "= h_" + i, "(mod p)")
        print(g_i, "^ y_" + i, "=", h_i, "( mod", p, ")")
        y_i = group(p, g_i, h_i, ord)
        # если решения нет, то возвращаем -1
        if y_i == -1:
            print("No solution")
            return -1
        y.append(y_i)
        print("y_" + i, "=", y_i, "( mod", ord, ")")
        i = int(i)
        i += 1
    # восстановление x
    # решаем систему сравнений с помощью КТО
    x = 0
    for i in range(len(y)):
        p_i, e_i = sorted(factors.items())[i]
        ord = pow(p_i, e_i)
        # N_i = N // ord
        N_i = gmpy2.divexact(mpz(N), mpz(ord))
        # x += l[i] * N_i * pow(N_i, -1, ord)
        N_i_inv = gmpy2.invert(mpz(N_i), mpz(ord))
        x = gmpy2.add(mpz(x), mpz(gmpy2.mul(mpz(y[i]), mpz(gmpy2.mul(mpz(N_i), mpz(N_i_inv))))))
    return x % N
