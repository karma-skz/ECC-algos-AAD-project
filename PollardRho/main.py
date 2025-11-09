"""
Pollard-rho for ECDLP (Floyd cycle detection, m-way partitions, negation map)
Reads 5-line positional input (p | a b | Gx Gy | n | Qx Qy) from a file path given as CLI arg,
or defaults to ./input/filename.txt relative to this file.
Outputs metadata, stats, and verification, aligned with other modules.
"""

from typing import Optional, Tuple, List
import random
import time
from pathlib import Path
import sys

Point = Optional[Tuple[int, int]]

# ---------- EC helpers ---------- #

def egcd(a: int, b: int) -> Tuple[int, int, int]:
    """Extended gcd: returns (g, x, y) such that a*x + b*y = g = gcd(a,b)."""
    if b == 0:
        return (abs(a), 1 if a >= 0 else -1, 0)
    else:
        g, x1, y1 = egcd(b, a % b)
        return (g, y1, x1 - (a // b) * y1)

def mod_inv_euclid(x: int, m: int) -> int:
    """Modular inverse using extended Euclid. Raises ZeroDivisionError if inverse doesn't exist."""
    x %= m
    g, u, _ = egcd(x, m)
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
    if P is None:
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

# ---------- Sanity checks ---------- #

def is_curve_valid(p: int, a: int, b: int) -> bool:
    # discriminant Δ = -16(4a^3 + 27b^2) ≠ 0  <=>  4a^3 + 27b^2 ≠ 0 mod p
    return ((4 * pow(a % p, 3, p)) + (27 * pow(b % p, 2, p))) % p != 0

def is_point_on_curve(P: Point, p: int, a: int, b: int) -> bool:
    if P is None:
        return True
    x, y = P
    if not (0 <= x < p and 0 <= y < p):
        return False
    return (y * y - (x * x * x + a * x + b)) % p == 0

# ---------- Rho helpers (negation map + partitions) ---------- #

def canonicalize(X: Point, A: int, B: int, p: int, n: int):
    """
    Apply the negation map to choose a canonical representative between X and -X.
    We flip to the representative with smaller y (ties broken implicitly).
    Adjust (A,B) accordingly: -(A*G + B*Q) = (-A)*G + (-B)*Q.
    """
    if X is None:
        return X, A, B
    x, y = X
    # choose the representative with y <= p - y
    if y > (p - y) % p:
        X = (x % p, (-y) % p)
        A = (-A) % n
        B = (-B) % n
    return X, A, B

def make_partition_table(a: int, p: int, G: Point, Q: Point, n: int, m: int = 16):
    """
    Build an m-way partition jump table: for each i, choose random (u_i, v_i),
    set R_i = u_i*G + v_i*Q. On bucket i, we do X <- X + R_i and (A,B) <- (A+u_i, B+v_i).
    """
    UV: List[Tuple[int, int]] = []
    Rtab: List[Point] = []
    for _ in range(m):
        while True:
            u = random.randrange(n)
            v = random.randrange(n)
            if (u | v) != 0:  # not both zero
                break
        Ri = ec_point_add(ec_scalar_mul(u, G, a, p), ec_scalar_mul(v, Q, a, p), a, p)
        UV.append((u, v))
        Rtab.append(Ri)
    return UV, Rtab

# ---------- Pollard-rho implementation ---------- #

def pollard_rho_ecdlp(p: int, a: int, b: int,
                      G: Point, Q: Point, n: int,
                      max_steps: int = 2_000_000,
                      partition_m: int = 16) -> Tuple[Optional[int], int]:
    """
    Pollard-rho for ECDLP using Floyd's cycle detection, m-way partitions, and negation map.
    Returns (d, steps) where d satisfies Q = d*G (mod n), else (None, steps).
    """
    if G is None or Q is None:
        return None, 0

    UV, Rtab = make_partition_table(a, p, G, Q, n, partition_m)

    def random_state():
        A = random.randrange(n)
        B = random.randrange(n)
        X = ec_point_add(ec_scalar_mul(A, G, a, p), ec_scalar_mul(B, Q, a, p), a, p)
        return canonicalize(X, A % n, B % n, p, n)

    def step(state):
        X, A, B = state
        idx = 0 if X is None else (X[0] % partition_m)
        Ri = Rtab[idx]
        ui, vi = UV[idx]
        Xn = ec_point_add(X, Ri, a, p)
        An = (A + ui) % n
        Bn = (B + vi) % n
        return canonicalize(Xn, An, Bn, p, n)

    # Floyd cycle detection
    tort = random_state()
    hare = step(tort)
    steps = 0

    while steps < max_steps:
        tort = step(tort)
        hare = step(step(hare))
        steps += 1

        Xt, At, Bt = tort
        Xh, Ah, Bh = hare

        if Xt == Xh:
            # We have (At*G + Bt*Q) == (Ah*G + Bh*Q)
            # -> (At - Ah)*G == (Bh - Bt)*Q == (Bh - Bt)*d*G
            # -> (At - Ah) == d*(Bh - Bt) (mod n)
            num = (At - Ah) % n
            den = (Bh - Bt) % n
            # Solve num = d * den (mod n)
            from math import gcd as math_gcd
            g = math_gcd(den, n)
            if g == 1:
                try:
                    den_inv = mod_inv_euclid(den, n)
                except ZeroDivisionError:
                    return None, steps
                d = (num * den_inv) % n
                if ec_scalar_mul(d, G, a, p) == Q:
                    return d, steps
                return None, steps
            else:
                # If consistent, solve modulo n/g and lift to g candidates.
                if num % g != 0:
                    return None, steps  # inconsistent; restart
                n_ = n // g
                den_ = (den // g) % n_
                num_ = (num // g) % n_
                try:
                    den_inv_ = mod_inv_euclid(den_, n_)
                except ZeroDivisionError:
                    return None, steps
                d0 = (num_ * den_inv_) % n_
                # Lift: d = d0 + k*(n/g), k=0..g-1
                for k in range(g):
                    d_candidate = (d0 + k * n_) % n
                    if ec_scalar_mul(d_candidate, G, a, p) == Q:
                        return d_candidate, steps
                return None, steps

    return None, steps

# ---------- Input loader ---------- #

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
    # Optional: make runs reproducible by setting a seed via env/arg if you like
    # random.seed(42)

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
        print("Warning: n*G != O; provided n may not be the exact order (use the prime order subgroup).")

    # Config
    method = "Pollard-rho (m-way partitions + negation map)"
    partition_m = 32       # try 16–64; 32 is a good default
    max_steps = 5_000_000  # per attempt
    attempts = 50          # independent restarts

    steps_total = 0
    start = time.perf_counter()
    d_found = None

    for attempt in range(1, attempts + 1):
        d_try, steps = pollard_rho_ecdlp(
            p, a, b, G, Q, n,
            max_steps=max_steps,
            partition_m=partition_m
        )
        steps_total += steps
        if d_try is not None:
            d_found = d_try
            break

    elapsed = time.perf_counter() - start

    # Output metadata & stats
    print(f"Curve: p={p}, a={a}, b={b}")
    print(f"G=({G[0]},{G[1]}), Q=({Q[0]},{Q[1]}), n={n}")
    print(f"Method: {method}")
    print(f"Partitions (m): {partition_m}")
    print(f"Steps taken (total): {steps_total}")
    print(f"Attempts: {attempts}")
    print(f"Time elapsed: {elapsed:.6f} s")

    if d_found is not None:
        print(f"Found d: {d_found}")
        Q_check = ec_scalar_mul(d_found % n, G, a, p)
        match = (Q_check == Q)
        print(f"Verified: {match}")
        print(f"Match: {match}")
    else:
        print("No solution found (None).")
        print("Tips: ensure n is the prime order of G; consider increasing m, attempts, or max_steps.")
