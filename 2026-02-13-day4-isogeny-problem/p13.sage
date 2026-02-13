import time
load("func.sage")

p = 14069954659
Fp2.<i> = GF(p**2, modulus=x^2 + 1)
j1 = random_ss_j_from_j0(Fp2(1728), ceil(log(p, 2)))
j2 = random_ss_j_from_j0(Fp2(1728), ceil(log(p, 2)))
print(f"j1 = {j1}")
print(f"j2 = {j2}")
print("Starting MITM...")
t0 = time.time()
chain_mitm = mitm_search(j1, j2, ceil(log(p, 2)))
t1 = time.time()
print(f"Time for MITM: {t1 - t0} seconds")