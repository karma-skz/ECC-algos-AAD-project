import os, random
from typing import Optional, Tuple

Point = Optional[Tuple[int, int]]

# -------- prime utilities -------- #

def mr_is_prime(n, k=12):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    d = n - 1
    r = 0
    while d % 2 == 0:
        d >>= 1
        r += 1
    for _ in range(k):
        a = random.randrange(2, n - 2)
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for __ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def rand_prime_40bit(d):
    """Generate a 40-bit prime p ≡ 11 mod 12 so Tonelli is fast."""
    target = 11  # p ≡ 11 (mod 12) => p%3=2 AND p%4=3
    while True:
        p = random.getrandbits(40) | 1
        if p <= d:
            continue
        if p % 12 != target:
            continue
        if p.bit_length() == 40 and mr_is_prime(p):
            return p

def egcd(a, b):
    if b == 0: return (a, 1, 0)
    g, x, y = egcd(b, a % b)
    return (g, y, x - (a // b) * y)

def mod_inv(a, n):
    g, x, _ = egcd(a % n, n)
    if g != 1:
        raise Exception("no inverse")
    return x % n

# -------- EC ops -------- #

def tonelli(n, p):
    """Fast path because p ≡ 3 (mod 4)."""
    if pow(n, (p - 1) // 2, p) != 1:
        return None
    return pow(n, (p + 1) // 4, p)

def ec_add(P, Q, a, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P != Q:
        lam = ((y2 - y1) * mod_inv(x2 - x1, p)) % p
    else:
        lam = ((3 * x1 * x1 + a) * mod_inv(2 * y1, p)) % p
    x3 = (lam * lam - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)

def ec_mul(k, P, a, p):
    R = None
    Q = P
    while k:
        if k & 1:
            R = ec_add(R, Q, a, p)
        Q = ec_add(Q, Q, a, p)
        k >>= 1
    return R

# -------- full-order test -------- #

def trial_factor(n):
    i = 2
    fac = {}
    N = n
    while i * i <= N:
        while N % i == 0:
            fac[i] = fac.get(i, 0) + 1
            N //= i
        i += 1 if i == 2 else 2
    if N > 1:
        fac[N] = fac.get(N, 0) + 1
    return fac

def has_full_order(G, n, a, p):
    fac = trial_factor(n)
    for q in fac.keys():
        if ec_mul(n // q, G, a, p) is None:
            return False
    return True

# -------- generator -------- #

def generate():
    os.makedirs("input", exist_ok=True)
    # d_secret = 9869284574
    d_secret = 696969

    for idx in range(1, 6):  # exactly 5 cases
        p = rand_prime_40bit(d_secret)
        a = 0
        b = 1
        n = p + 1  # order of curve for p ≡ 2 mod 3

        # find a full-order point
        while True:
            x = random.randrange(0, p)
            rhs = (x * x * x + 1) % p
            y = tonelli(rhs, p)
            if y is None:
                continue
            G = (x, y)
            if has_full_order(G, n, a, p):
                break

        d_mod = d_secret % n
        Q = ec_mul(d_mod, G, a, p)

        # write files
        with open(f"input/testcase_{idx}.txt", "w") as f:
            f.write(f"{p}\n")
            f.write(f"{a} {b}\n")
            f.write(f"{G[0]} {G[1]}\n")
            f.write(f"{n}\n")
            f.write(f"{Q[0]} {Q[1]}\n")

        with open(f"input/answer_{idx}.txt", "w") as f:
            f.write(str(d_mod))

        print(f"[✓] testcase_{idx} generated (p={p}, n={n}, d={d_mod})")

if __name__ == "__main__":
    generate()