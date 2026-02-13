load("func.sage")

p = 161531
Fp2.<i> = GF(p**2, modulus=x^2 + 1)

j1 = random_ss_j_from_j0(Fp2(1728), 100)
j2 = random_ss_j_from_j0(Fp2(1728), 100)

path, degs = GD_search(j1, j2, 30, [3, 5, 7], 5)
assert path[0] == j1
assert path[-1] == j2
for i in range(len(degs)):
    assert classical_modular_polynomial(degs[i])(path[i], path[i + 1]) == 0
print("GD search succeeded")
print(f"Path length: {len(path) - 1}")
print(f"path: {path}")
print(f"degs: {degs}")