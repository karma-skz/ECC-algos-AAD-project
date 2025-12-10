"""
Pohlig-Hellman with Residue Leakage - BONUS
Adaptation: Skips computation for leaked sub-moduli.
"""
import sys, time, ctypes, math
from collections import defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import EllipticCurve, Point, load_input, crt_combine
from utils.bonus_utils import print_bonus_result

USE_CPP = False
ecc_lib = None
lib_paths = [Path(__file__).parent.parent / "utils" / "cpp" / "ecc_fast.so"]
for p in lib_paths:
    if p.exists():
        ecc_lib = ctypes.CDLL(str(p))
        ecc_lib.scalar_mult.argtypes = [ctypes.c_longlong]*6 + [ctypes.POINTER(ctypes.c_longlong)]*2
        ecc_lib.scalar_mult.restype = ctypes.c_int
        USE_CPP = True
        break

def fast_mult(k, G, curve):
    if not USE_CPP: return curve.scalar_multiply(k, G)
    rx, ry = ctypes.c_longlong(), ctypes.c_longlong()
    # FIX: No % p
    valid = ecc_lib.scalar_mult(k, G[0], G[1], curve.a, curve.b, curve.p, ctypes.byref(rx), ctypes.byref(ry)) #type:ignore
    return (rx.value, ry.value) if valid else None

def trial_factor(n):
    """Simple factorization using trial division."""
    factors = defaultdict(int)
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] += 1
            n //= d
        d += 1
    if n > 1: factors[n] += 1
    return dict(factors)

def bsgs_small(curve, G, Q, order):
    """Baby-step giant-step for small order subgroups using C++ optimization."""
    if Q is None:
        return 0
    
    m = int(math.isqrt(order)) + 1
    
    # Baby steps: Store j*G for j = 0, 1, ..., m-1
    baby_table = {}
    R = None
    
    for j in range(m):
        if R not in baby_table:
            baby_table[R] = j
        R = curve.add(R, G)
    
    # Giant steps
    mG = fast_mult(m, G, curve)
    neg_mG = curve.negate(mG) if mG else None
    
    gamma = Q
    for i in range(m):
        if gamma in baby_table:
            return i * m + baby_table[gamma]
        gamma = curve.add(gamma, neg_mG)
    
    return None

def solve_with_leak(curve, G, Q, n, leaks={}):
    factors = trial_factor(n)
    congruences = []
    
    print(f"{'MODULUS':<15} | {'SOURCE':<15} | {'TIME'}")
    print(f"{'-'*45}")
    
    for q, e in factors.items():
        mod = q**e
        t0 = time.perf_counter()
        
        if mod in leaks:
            d_i = leaks[mod]
            src = "LEAKED"
        else:
            h = n // mod
            d_i = bsgs_small(curve, fast_mult(h, G, curve), fast_mult(h, Q, curve), mod)
            src = "Computed"
            
        print(f"{mod:<15} | {src:<15} | {time.perf_counter()-t0:.6f}s")
        if d_i is not None: congruences.append((d_i, mod))
        
    d, _ = crt_combine(congruences)
    return d % n

def main():
    if len(sys.argv) < 2: return
    p, a, b, G, n, Q = load_input(Path(sys.argv[1]))
    curve = EllipticCurve(a, b, p)
    try:
        num = Path(sys.argv[1]).stem.split('_')[1]
        with open(Path(sys.argv[1]).parent / f"answer_{num}.txt") as f: d_real = int(f.read())
    except: d_real = 12345

    print(f"\n{'='*60}")
    print(f"POHLIG-HELLMAN LEAKAGE DEMO")
    print(f"{'='*60}")
    
    factors = trial_factor(n)
    if not factors: return
    
    # Leak largest factor
    q, e = sorted(factors.items())[-1]
    mod = q**e
    leaks = {mod: d_real % mod}
    
    print(f"Simulating leak of d mod {mod}")
    t0 = time.perf_counter()
    d = solve_with_leak(curve, G, Q, n, leaks)
    elapsed = time.perf_counter() - t0
    print(f"{'='*60}")
    print(f"Result: {'[SUCCESS]' if d==d_real else '[FAILED]'}")

    print_bonus_result("PohligHellman", "success" if d == d_real else "fail", elapsed, 0, {"leaked_modulus": mod})

if __name__ == "__main__":
    main()