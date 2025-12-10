"""
Optimized Las Vegas (Birthday Attack) with C++ backend.

Path Correction:
    Loads C++ library from: ../utils/cpp/ecc_fast.so
    
Feature:
    Includes a "Self-Destruct Timer" to prevent hanging the comparison script.
"""

import sys
import time
import random
import ctypes
import math
from pathlib import Path
from typing import Optional, Tuple, Dict

# Add parent directory to path so we can import utils
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import EllipticCurve, Point, load_input, mod_inv

# -----------------------------------------------------------------------------
# C++ Library Loading
# -----------------------------------------------------------------------------
USE_CPP = False
ecc_lib = None
lib_path = Path(__file__).parent.parent / "utils" / "cpp" / "ecc_fast.so"

if lib_path.exists():
    try:
        ecc_lib = ctypes.CDLL(str(lib_path))
        ecc_lib.scalar_mult.argtypes = [
            ctypes.c_longlong, ctypes.c_longlong, ctypes.c_longlong,
            ctypes.c_longlong, ctypes.c_longlong, ctypes.c_longlong,
            ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_longlong)
        ]
        ecc_lib.scalar_mult.restype = ctypes.c_int
        USE_CPP = True
        print(f"✓ Using C++ optimized scalar multiplication")
    except Exception as e:
        print(f"⚠ Found {lib_path.name} but failed to load: {e}")
        USE_CPP = False
else:
    # Silent fallback to avoid spamming output in comparison script
    # print(f"⚠ C++ library not found at: {lib_path}")
    USE_CPP = False


def fast_scalar_mult(k: int, G: Point, curve: EllipticCurve) -> Optional[Point]:
    if not USE_CPP or G is None:
        return curve.scalar_multiply(k, G)
    result_x = ctypes.c_longlong()
    result_y = ctypes.c_longlong()
    is_valid = ecc_lib.scalar_mult( #type: ignore
        k % curve.p, G[0], G[1], curve.a, curve.b, curve.p,
        ctypes.byref(result_x), ctypes.byref(result_y)
    )
    return (result_x.value, result_y.value) if is_valid else None

def fast_point_add(P: Point, Q: Point, curve: EllipticCurve) -> Point:
    return curve.add(P, Q)


# -----------------------------------------------------------------------------
# ALGORITHM: Birthday Attack (Las Vegas)
# -----------------------------------------------------------------------------

def las_vegas_optimized(curve: EllipticCurve, G: Point, Q: Point, n: int, 
                       max_attempts: int = 5000000, 
                       timeout_sec: float = 20.0) -> Tuple[Optional[int], int]:
    """
    Optimized Las Vegas with Time Limit.
    Exits early if execution exceeds timeout_sec.
    """
    visited: Dict[Point, Tuple[int, int]] = {}
    
    # Dynamic Limit:
    # 1. Calculate Birthday Bound: sqrt(pi*n/2)
    # 2. Multiply by 10 to ensure >99.9% probability
    # 3. Cap at max_attempts (RAM safety)
    expected_points = int(math.sqrt(math.pi * n / 2))
    limit = min(max_attempts, max(2000, expected_points * 10))
    
    if n > 10000:
        print(f"  Target: Collecting ~{expected_points} points...", flush=True)

    start_time = time.time()

    for i in range(1, limit + 1):
        # 1. Self-Destruct Check (Every 1000 steps)
        if i % 1000 == 0:
            if (time.time() - start_time) > timeout_sec:
                # print(f"  [TIMEOUT] Exceeded {timeout_sec}s limit.")
                return None, i
            if i % 50000 == 0:
                print(f"  stored {i} points...", flush=True)

        # 2. Pick random coefficients
        a = random.randrange(n)
        b = random.randrange(n)
        
        # 3. Compute R = aG + bQ
        term1 = fast_scalar_mult(a, G, curve)
        term2 = fast_scalar_mult(b, Q, curve)
        R = fast_point_add(term1, term2, curve)
        
        if R is None: continue 
            
        # 4. Check for Collision
        if R in visited:
            a_old, b_old = visited[R]
            num = (a - a_old) % n
            den = (b_old - b) % n
            if den == 0: continue 

            g = math.gcd(den, n)
            
            # Case 1: Coprime (Standard)
            if g == 1:
                try:
                    d = (num * mod_inv(den, n)) % n
                    if fast_scalar_mult(d, G, curve) == Q: return d, i
                except: pass
            
            # Case 2: Composite Order (Subgroup Attack)
            elif num % g == 0:
                try:
                    n_prime = n // g
                    d_base = ((num // g) * mod_inv(den // g, n_prime)) % n_prime
                    for k in range(g):
                        d_cand = d_base + k * n_prime
                        if fast_scalar_mult(d_cand, G, curve) == Q: return d_cand, i
                except: pass
        
        # 5. Store point
        visited[R] = (a, b)
        
    return None, limit


def main():
    script_dir = Path(__file__).parent
    
    if len(sys.argv) > 1:
        input_path = Path(sys.argv[1])
    else:
        input_path = script_dir.parent / 'input' / 'testcase_1.txt'
    
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)
    
    try:
        p, a, b, G, n, Q = load_input(input_path)
        curve = EllipticCurve(a, b, p)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Solving ECDLP using Optimized Las Vegas (Birthday Attack)...")
    print(f"Curve: y^2 = x^3 + {a}x + {b} (mod {p})")
    print(f"G = ({G[0]}, {G[1]}), Q = ({Q[0]}, {Q[1]}), n = {n}") #type: ignore
    
    if n > 2**30:
        print("\n[WARNING] Curve > 30 bits. RAM usage will be high.")
    
    # -----------------------------------------------------
    # SETTINGS:
    # -----------------------------------------------------
    max_steps = 5_000_000 
    
    # TIMEOUT: Exit after 20s so run_comparison doesn't hang
    TIMEOUT_LIMIT = 20.0 
    
    start_time = time.perf_counter()
    
    d, attempts = las_vegas_optimized(curve, G, Q, n, max_steps, timeout_sec=TIMEOUT_LIMIT)
    
    elapsed = time.perf_counter() - start_time
    
    if d is not None:
        print(f"\n{'='*50}")
        print(f"Solution: d = {d}")
        print(f"Points Stored: {attempts}")
        print(f"Time: {elapsed:.6f} seconds")
        
        Q_verify = fast_scalar_mult(d, G, curve)
        verified = (Q_verify == Q)
        print(f"Verification: {'PASSED' if verified else 'FAILED'}")
        print(f"{'='*50}")
        if not verified: sys.exit(1)
    else:
        print(f"\nNo solution found after storing {attempts} points.")
        print(f"Time: {elapsed:.6f} seconds")
        sys.exit(1)

if __name__ == "__main__":
    main()