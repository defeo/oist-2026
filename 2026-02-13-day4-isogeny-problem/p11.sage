load("func.sage")

for p in Primes():
    if p < 100:
        continue
    if p > 200:
        break
    Fp2 = GF(p**2)
    n = num_of_codomain_j(Fp2, ceil(log(p, 2)))
    nss = floor(p/12)
    if p % 3 == 2:
        nss += 1
    if p % 4 == 3:
        nss += 1
    print(f"p = {p}, {n == nss}")