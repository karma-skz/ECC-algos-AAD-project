"""Pohlig-Hellman ECDLP - Optimized with C++"""
import sys, time, ctypes, math
from pathlib import Path
from collections import defaultdict
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import EllipticCurve, Point, load_input, mod_inv

# Load C++ library
USE_CPP = False
try:
    lib_path = Path(__file__).parent.parent / 'utils' / 'cpp' / 'ecc_fast.so'
    if lib_path.exists():
        ecc_lib = ctypes.CDLL(str(lib_path))
        ecc_lib.scalar_mult.argtypes = [ctypes.c_int64]*6 + [ctypes.POINTER(ctypes.c_int64)]*2
        ecc_lib.scalar_mult.restype = ctypes.c_int
        USE_CPP = True
        print("✓ Using C++ optimization")
except: pass

def fast_scalar_mult(k, G, curve):
    if not USE_CPP or G is None: return curve.scalar_multiply(k, G)
    rx, ry = ctypes.c_int64(), ctypes.c_int64()
    valid = ecc_lib.scalar_mult(ctypes.c_int64(k), ctypes.c_int64(G[0]), ctypes.c_int64(G[1]),
                                  ctypes.c_int64(curve.a), ctypes.c_int64(curve.b), ctypes.c_int64(curve.p),
                                  ctypes.byref(rx), ctypes.byref(ry))
    return None if valid == 0 else (rx.value, ry.value)

def factor(n):
    """Simple factorization."""
    factors = defaultdict(int)
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] += 1
            n //= d
        d += 1
    if n > 1: factors[n] += 1
    return dict(factors)

def bsgs_subproblem(curve, G_i, Q_i, p_i):
    """Baby-step giant-step for subproblem."""
    if G_i is None or Q_i is None:
        return None
    
    m = int(math.ceil(math.sqrt(p_i)))
    
    # Baby steps: store j*G_i
    baby_table = {}
    baby_table[None] = 0  # Identity point
    R = G_i
    for j in range(1, m + 1):
        if R is None:
            break
        baby_table[R] = j
        next_R = curve.add(R, G_i)
        if next_R == R:  # Point has small order
            break
        R = next_R
    
    # Giant steps
    mG = fast_scalar_mult(m, G_i, curve)
    if mG is None:
        # Try smaller search
        for j in range(p_i):
            if fast_scalar_mult(j, G_i, curve) == Q_i:
                return j
        return None
    
    neg_mG = (mG[0], (-mG[1]) % curve.p)
    
    gamma = Q_i
    for i in range(m + 1):
        if gamma in baby_table:
            return (i * m + baby_table[gamma]) % p_i
        if gamma is None:
            if None in baby_table:
                return (i * m + baby_table[None]) % p_i
        gamma = curve.add(gamma, neg_mG) if gamma else neg_mG
    
    return None

def crt(residues, moduli):
    """Chinese Remainder Theorem."""
    total = 0
    prod = 1
    for m in moduli: prod *= m
    for r, m in zip(residues, moduli):
        p = prod // m
        total += r * mod_inv(p, m) * p
    return total % prod

if __name__ == "__main__":
    p, a, b, G, n, Q = load_input(Path(sys.argv[1]) if len(sys.argv) > 1 else Path('input/test_1.txt'))
    curve = EllipticCurve(a, b, p)
    print(f"Pohlig-Hellman ECDLP: p={p}, n={n}")
    
    start = time.time()
    factors = factor(n)
    print(f"Factorization: {factors}")
    
    residues, moduli = [], []
    for prime, exp in factors.items():
        p_e = prime ** exp
        G_i = fast_scalar_mult(n // p_e, G, curve)
        Q_i = fast_scalar_mult(n // p_e, Q, curve)
        print(f"Solving for p^e = {prime}^{exp}...")
        d_i = bsgs_subproblem(curve, G_i, Q_i, p_e)
        if d_i is None:
            print(f"✗ Subproblem failed (order might be incorrect)")
            # Try brute force as fallback for small primes
            if p_e <= 1000:
                for attempt in range(p_e):
                    if fast_scalar_mult(attempt, G_i, curve) == Q_i:
                        d_i = attempt
                        print(f"  → Found via brute force: {d_i}")
                        break
            if d_i is None:
                print("✗ PohligHellman requires exact point order (not available)")
                print("  Note: n=p+1 is approximation, actual order may differ")
                sys.exit(1)
        residues.append(d_i)
        moduli.append(p_e)
    
    d = crt(residues, moduli)
    elapsed = time.time() - start
    
    if fast_scalar_mult(d, G, curve) == Q:
        print(f"✓ Solution: d = {d}")
        print(f"Time: {elapsed:.6f}s")
        print("Verification: PASSED")
        sys.exit(0)
    print("✗ Verification failed")
    sys.exit(1)
