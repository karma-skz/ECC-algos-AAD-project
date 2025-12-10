#!/usr/bin/env python3
"""
ECC ECDLP Test Case Generator (Prime-Subgroup Mode) - single file

- Finds primes p with p ≡ 11 (mod 12)
- Ensures q = (p + 1) / 12 is prime
- Constructs curve y^2 = x^3 + a*x + b (default a = 0)
- Finds generator G of exact prime order q via cofactor clearing + verification
- Writes test_cases/<bits>bit/case_i.txt and answer_i.txt
"""
import sys
import random
from pathlib import Path
import traceback

# --- Configuration ---
MIN_BITS = 10
MAX_BITS = 60

# -------------------------
# Elliptic curve (short Weierstrass)
# -------------------------
class EllipticCurve:
    def __init__(self, a, b, p):
        self.a = a % p
        self.b = b % p
        self.p = p

    def add(self, P, Q):
        # None is point at infinity
        if P is None:
            return Q
        if Q is None:
            return P

        x1, y1 = P
        x2, y2 = Q
        p = self.p

        # P == -Q -> infinity
        if x1 == x2 and (y1 + y2) % p == 0:
            return None

        if P != Q:
            den = (x2 - x1) % p
            if den == 0:
                return None
            lam = ((y2 - y1) * pow(den, -1, p)) % p
        else:
            # Doubling; if y1 == 0 -> tangent vertical -> infinity
            if y1 % p == 0:
                return None
            den = (2 * y1) % p
            lam = ((3 * x1 * x1 + self.a) * pow(den, -1, p)) % p

        x3 = (lam * lam - x1 - x2) % p
        y3 = (lam * (x1 - x3) - y1) % p
        return (x3, y3)

    def scalar_multiply(self, k, P):
        R = None
        Q = P
        while k > 0:
            if k & 1:
                R = self.add(R, Q)
            Q = self.add(Q, Q)
            k >>= 1
        return R

# -------------------------
# Primality test (Miller-Rabin)
# -------------------------
def is_prime_miller_rabin(n, k=20):
    if n in (2, 3):
        return True
    if n <= 1 or n % 2 == 0:
        return False
    # write n-1 = 2^r * d
    r = 0
    d = n - 1
    while d % 2 == 0:
        r += 1
        d //= 2
    for _ in range(k):
        a = random.randrange(2, n - 1)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

# -------------------------
# Find p (p ≡ 11 mod 12) and q = (p+1)/12 prime
# -------------------------
def find_strong_prime_pair(bits, seed_offset=0, max_attempts=200000):
    """
    Find prime p with p ≡ 11 (mod 12) and q = (p+1)/12 prime.
    Uses deterministic search with random offset for diversity.
    Returns (p, q) or (None, None).
    """
    random.seed(bits * 10000 + seed_offset)
    min_val = 1 << (bits - 1)
    max_val = (1 << bits) - 1
    
    # For small bit sizes, use a more deterministic approach
    # Start at a random offset within the range but search more systematically
    range_size = max_val - min_val
    offset = random.randint(0, range_size) if range_size > 12 else 0
    
    # Calculate starting candidate
    start = min_val + offset
    candidate = start - (start % 12) + 11
    if candidate < min_val:
        candidate += 12
    
    # First pass: search forward from starting point
    temp = candidate
    attempts = 0
    while attempts < max_attempts // 2 and temp <= max_val:
        if is_prime_miller_rabin(temp):
            q = (temp + 1) // 12
            if q > 1 and is_prime_miller_rabin(q):
                return temp, q
        temp += 12
        attempts += 1
    
    # Second pass: search backward from starting point (wrap around)
    temp = candidate - 12
    while attempts < max_attempts and temp >= min_val:
        if is_prime_miller_rabin(temp):
            q = (temp + 1) // 12
            if q > 1 and is_prime_miller_rabin(q):
                return temp, q
        temp -= 12
        attempts += 1
    
    # If still not found, do exhaustive search from beginning
    temp = min_val - (min_val % 12) + 11
    if temp < min_val:
        temp += 12
    while temp <= max_val:
        if is_prime_miller_rabin(temp):
            q = (temp + 1) // 12
            if q > 1 and is_prime_miller_rabin(q):
                return temp, q
        temp += 12
    
    return None, None

# -------------------------
# Quadratic residue check and sqrt (p % 4 == 3 branch)
# -------------------------
def is_quadratic_residue(rhs, p):
    v = pow(rhs % p, (p - 1) // 2, p)
    return v == 1 or v == 0

def sqrt_mod_p(rhs, p):
    rhs = rhs % p
    if rhs == 0:
        return 0
    if p % 4 != 3:
        raise ValueError("sqrt_mod_p only implemented for p % 4 == 3")
    return pow(rhs, (p + 1) // 4, p)

# -------------------------
# Find subgroup generator of exact order q
# -------------------------
def find_subgroup_generator(curve, p, q, cofactor=12, max_tries=20000):
    """
    Returns G such that order(G) == q (checks q*G == O and G != O).
    Raises RuntimeError if not found.
    """
    if p % 4 != 3:
        raise RuntimeError("p % 4 != 3; sqrt_mod_p not available in this script")

    for _ in range(max_tries):
        x = random.randint(0, p - 1)
        rhs = (pow(x, 3, p) + curve.a * x + curve.b) % p

        if not is_quadratic_residue(rhs, p):
            continue

        try:
            y = sqrt_mod_p(rhs, p)
        except ValueError:
            continue

        P = (x, y)
        G = curve.scalar_multiply(cofactor, P)
        if G is None:
            continue

        # verify exact order: q * G == O (None) and G != O
        if curve.scalar_multiply(q, G) is None:
            return G

    raise RuntimeError("Failed to find subgroup generator after many tries")

# -------------------------
# Generate test cases
# -------------------------
def generate_test_cases_for_bits(k, num_cases=5):
    if k < MIN_BITS or k > MAX_BITS:
        raise ValueError(f"Bit length {k} out of range ({MIN_BITS}-{MAX_BITS})")

    created = 0
    for case_num in range(1, num_cases + 1):
        p, q = find_strong_prime_pair(k, seed_offset=case_num)
        if not p:
            print(f"  [Skip] case {case_num}: no p,q found for {k}-bit")
            continue

        a = 0
        b = random.randint(1, p - 1)
        curve = EllipticCurve(a, b, p)

        try:
            G = find_subgroup_generator(curve, p, q, cofactor=12)
        except Exception as e:
            print(f"  [Skip] case {case_num}: generator search failed: {e}")
            continue

        d = random.randint(1, q - 1) #type:ignore
        Q = curve.scalar_multiply(d, G)
        if Q is None:
            print(f"  [Skip] case {case_num}: Q == O (rare), retry")
            continue

        test_dir = Path(__file__).parent / 'test_cases' / f'{k:02d}bit'
        test_dir.mkdir(parents=True, exist_ok=True)

        with open(test_dir / f'case_{case_num}.txt', 'w') as f:
            f.write(f'{p}\n{a} {b}\n{G[0]} {G[1]}\n{q}\n{Q[0]} {Q[1]}\n')

        with open(test_dir / f'answer_{case_num}.txt', 'w') as f:
            f.write(f'{d}\n')

        print(f"  ✓ case {case_num}: p={p} (q bits ~= {q.bit_length()})") #type:ignore
        created += 1

    return created

# -------------------------
# Small validator utility (optional)
# -------------------------
def validate_case_file(case_path):
    """
    Quick check: parse case file and ensure q*G == O and Q = d*G using answer file.
    Returns True if passes, False otherwise.
    """
    data = case_path.read_text().strip().splitlines()
    if len(data) < 5:
        return False
    p = int(data[0].strip())
    a, b = map(int, data[1].split())
    Gx, Gy = map(int, data[2].split())
    q = int(data[3].strip())
    Qx, Qy = map(int, data[4].split())

    curve = EllipticCurve(a, b, p)
    G = (Gx, Gy)
    Q = (Qx, Qy)

    # check q*G == O
    if curve.scalar_multiply(q, G) is not None:
        return False

    # answer file
    ans_path = case_path.with_name('answer_' + case_path.name.split('_')[1])
    if not ans_path.exists():
        return False
    d = int(ans_path.read_text().strip())
    Qtest = curve.scalar_multiply(d, G)
    return Qtest == Q

# -------------------------
# CLI
# -------------------------
def main():
    args = sys.argv[1:]
    start_bit = 10
    end_bit = 20
    cases_per_bit = 5

    if len(args) >= 1:
        start_bit = int(args[0])
    if len(args) >= 2:
        end_bit = int(args[1])
    if len(args) >= 3:
        cases_per_bit = int(args[2])

    print("="*60)
    print("ECC Generator: Prime-Subgroup Mode (Fixed)")
    print("Generates curves with order N = 12*q where q is prime.")
    print("="*60)

    for k in range(start_bit, end_bit + 1):
        try:
            print(f"Generating for {k}-bit (attempting {cases_per_bit} cases)...")
            created = generate_test_cases_for_bits(k, cases_per_bit)
            print(f"  -> created {created} cases for {k}-bit")
        except Exception:
            print(f"Error generating cases for {k}-bit:")
            traceback.print_exc()

    print("\nDone.")

if __name__ == "__main__":
    main()
