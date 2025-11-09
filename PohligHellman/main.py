"""
Pohlig-Hellman for ECDLP
Reads 5-line positional input (p | a b | Gx Gy | n | Qx Qy) from a file path given as CLI arg,
or defaults to ./input/filename.txt relative to this file.
Outputs metadata, stats, and verification, aligned with other modules.
"""

from typing import Optional, Tuple, Dict, List
from pathlib import Path
import sys
import time
import math

Point = Optional[Tuple[int, int]]

# ------------------------- EC helpers ------------------------- #

def egcd(a: int, b: int) -> Tuple[int, int, int]:
    """Extended gcd: returns (g, x, y) with a*x + b*y = g."""
    if b == 0:
        return (a, 1, 0)
    g, x1, y1 = egcd(b, a % b)
    return (g, y1, x1 - (a // b) * y1)

def crt_combine(congs: List[Tuple[int, int]]) -> Tuple[int, int]:
    """
    Combine list of congruences (r_i mod m_i) using CRT.
    Returns (r, M) where r is solution modulo M = product(m_i).
    Assumes moduli are pairwise coprime.
    """
    r = 0
    M = 1
    for _, m in congs:
        M *= m
    for a_i, m_i in congs:
        Mi = M // m_i
        g, inv, _ = egcd(Mi, m_i)
        if g != 1:
            raise ValueError("Moduli not coprime in CRT combine")
        r = (r + a_i * Mi * (inv % m_i)) % M
    return r % M, M

def mod_inv_euclid(x: int, m: int) -> int:
    """Modular inverse using extended Euclid."""
    x %= m
    g, u, v = egcd(x, m)
    if g != 1:
        raise ZeroDivisionError(f"no inverse: gcd({x},{m}) = {g}")
    return u % m

def ec_point_neg(P: Point, p: int) -> Point:
    if P is None:
        return None
    x, y = P
    return (x % p, (-y) % p)

def ec_point_add(P: Point, Q: Point, a: int, p: int) -> Point:
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P != Q:
        num = (y2 - y1) % p
        den = (x2 - x1) % p
        lam = (num * mod_inv_euclid(den, p)) % p
    else:
        if y1 % p == 0:
            return None
        num = (3 * (x1 * x1 % p) + a) % p
        den = (2 * y1) % p
        lam = (num * mod_inv_euclid(den, p)) % p
    x3 = (lam * lam - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)

def ec_scalar_mul(k: int, P: Point, a: int, p: int) -> Point:
    """Double-and-add scalar multiplication."""
    # BUG FIX: previously checked k % p == 0 which is incorrect; should check k == 0.
    if P is None or k == 0:
        return None
    if k < 0:
        return ec_point_neg(ec_scalar_mul(-k, P, a, p), p)
    R: Point = None
    Q: Point = P
    while k:
        if k & 1:
            R = ec_point_add(R, Q, a, p)
        Q = ec_point_add(Q, Q, a, p)
        k >>= 1
    return R

# ------------------------- Utilities ------------------------- #

def is_curve_valid(p: int, a: int, b: int) -> bool:
    return ((4 * (a % p) ** 3) + (27 * (b % p) ** 2)) % p != 0

def is_point_on_curve(P: Point, p: int, a: int, b: int) -> bool:
    if P is None:
        return True
    x, y = P
    if not (0 <= x < p and 0 <= y < p):
        return False
    return (y * y - (x * x * x + a * x + b)) % p == 0

def trial_factor(n: int) -> Dict[int, int]:
    """Simple trial division factorization (sufficient for small n in toy examples)."""
    i = 2
    fac: Dict[int, int] = {}
    N = n
    while i * i <= N:
        while N % i == 0:
            fac[i] = fac.get(i, 0) + 1
            N //= i
        i += 1 if i == 2 else 2  # skip even numbers after 2
    if N > 1:
        fac[N] = fac.get(N, 0) + 1
    return fac

# ------------------------- Small BSGS (works for small orders q) ------------------------- #

def bsgs_discrete_log(p: int, a: int, G: Point, Q: Point, order: int) -> Optional[int]:
    """
    Baby-step Giant-step to find x in [0, order-1] with x*G = Q.
    order should be small (prime factor q) — used inside Pohlig-Hellman.
    """
    if Q is None:
        return 0
    m = int(math.isqrt(order)) + 1
    baby: Dict[Point, int] = {}
    R: Point = None
    # baby steps: store j such that R = j*G (R starts at O -> j=0)
    for j in range(m):
        if R not in baby:
            baby[R] = j
        R = ec_point_add(R, G, a, p)  # advance by G

    # compute -m*G
    mG = ec_scalar_mul(m, G, a, p)
    neg_mG = ec_point_neg(mG, p)

    Gamma = Q
    for i in range(m + 1):
        if Gamma in baby:
            j = baby[Gamma]
            x = (i * m + j) % order
            return x
        Gamma = ec_point_add(Gamma, neg_mG, a, p)
    return None

# ------------------------- Pohlig-Hellman ------------------------- #

def pohlig_hellman_ecdlp(p: int, a: int, b: int, G: Point, Q: Point, n: int) -> Tuple[Optional[int], Dict[int, int], int]:
    """
    Pohlig-Hellman solver.
    Returns (d, factorization, total_bsgs_calls).
    d is discrete log modulo n (0..n-1) if found, else None.
    factorization is dict of prime->exponent.
    total_bsgs_calls counts the number of small BSGS solves performed.
    """
    # factor n
    fac = trial_factor(n)
    if len(fac) == 0:
        fac = {n: 1}

    congruences: List[Tuple[int, int]] = []  # (residue, modulus)
    total_bsgs = 0

    # For each prime power q^e
    for q, e in fac.items():
        n_i = q ** e
        h = n // n_i
        # Lifted points so that G1 has order dividing n_i
        G1 = ec_scalar_mul(h, G, a, p)
        Q1 = ec_scalar_mul(h, Q, a, p)

        # compute d mod q^e
        x_mod = 0  # current recovered value modulo q^k
        for k in range(e):
            # t = q^{e-1-k}
            t = q ** (e - 1 - k)
            # tmp = Q1 - x_mod * G1
            tmp = ec_point_add(Q1, ec_point_neg(ec_scalar_mul(x_mod, G1, a, p), p), a, p)
            # c_k = t * tmp
            c_k = ec_scalar_mul(t, tmp, a, p)
            # g_k = t * G1  -> has order q
            g_k = ec_scalar_mul(t, G1, a, p)
            # solve d_k in [0, q-1] with d_k * g_k = c_k
            d_k = bsgs_discrete_log(p, a, g_k, c_k, q)
            total_bsgs += 1
            if d_k is None:
                return None, fac, total_bsgs
            # update x_mod += d_k * q^k
            x_mod = (x_mod + d_k * (q ** k)) % n_i

        congruences.append((x_mod, n_i))

    # Combine congruences with CRT (moduli should be coprime)
    try:
        d, M = crt_combine(congruences)
    except ValueError:
        return None, fac, total_bsgs

    # canonicalize d modulo n
    d = d % n
    return d, fac, total_bsgs

# ------------------------- Input loader ------------------------- #

def load_positional_input(file_path: Path):
    with file_path.open('r') as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    if len(lines) < 5:
        raise ValueError("Input file must contain 5 non-empty lines: p | a b | Gx Gy | n | Qx Qy")
    def ints(s: str):
        return list(map(int, s.split()))
    p = ints(lines[0])[0]
    a, b = ints(lines[1])
    Gx, Gy = ints(lines[2])
    n = ints(lines[3])[0]
    Qx, Qy = ints(lines[4])
    return p, a, b, (Gx, Gy), n, (Qx, Qy)

# ------------------------- Runner ------------------------- #
if __name__ == "__main__":
    script_dir = Path(__file__).parent
    default_path = script_dir / 'input' / 'filename.txt'
    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else default_path
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    p, a, b, G, n, Q = load_positional_input(input_path)

    # Sanity checks
    if not is_curve_valid(p, a, b):
        raise ValueError("Invalid curve: 4a^3 + 27b^2 ≡ 0 (mod p)")
    if G is None or Q is None:
        raise ValueError("G/Q cannot be point at infinity")
    if not is_point_on_curve(G, p, a, b):
        raise ValueError("G is not on the curve")
    if not is_point_on_curve(Q, p, a, b):
        raise ValueError("Q is not on the curve")
    nG = ec_scalar_mul(n, G, a, p)
    if nG is not None:
        print("Warning: n*G != O; provided n may not be the exact order")

    method = "Pohlig-Hellman"
    start = time.perf_counter()
    d_found, fac, calls = pohlig_hellman_ecdlp(p, a, b, G, Q, n)
    elapsed = time.perf_counter() - start

    # Output metadata & stats
    print(f"Curve: p={p}, a={a}, b={b}")
    print(f"G=({G[0]},{G[1]}), Q=({Q[0]},{Q[1]}), n={n}")
    print(f"Method: {method}")
    print(f"Factorization of n: {fac}")
    print(f"BSGS calls (small q solves): {calls}")
    print(f"Time elapsed: {elapsed:.6f} s")

    if d_found is not None:
        print(f"Found d: {d_found}")
        Q_check = ec_scalar_mul(d_found % n, G, a, p)
        match = (Q_check == Q)
        print(f"Verified: {match}")
        print(f"Match: {match}")
    else:
        print("No solution found (None)")
