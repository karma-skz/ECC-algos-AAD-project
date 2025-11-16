"""
Pollard's Rho with Partial Key Leakage - BONUS Implementation (C++ Optimized)

Demonstrates how bit leaks and intervals restrict the random walk.
Uses C++ backend for fast operations.
"""

import sys
import time
import ctypes
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import (EllipticCurve, Point, load_input, KeyLeakage,
                   candidates_from_lsb_leak, candidates_from_msb_leak,
                   format_leak_info, calculate_speedup, calculate_search_reduction)

# Load C++ library
USE_CPP = False
try:
    lib_path = Path(__file__).parent / 'ecc_fast.so'
    if lib_path.exists():
        ecc_lib = ctypes.CDLL(str(lib_path))
        ecc_lib.scalar_mult.argtypes = [ctypes.c_int64]*6 + [ctypes.POINTER(ctypes.c_int64)]*2
        ecc_lib.scalar_mult.restype = ctypes.c_int
        ecc_lib.point_add.argtypes = [ctypes.c_int64]*8 + [ctypes.POINTER(ctypes.c_int64)]*2
        ecc_lib.point_add.restype = ctypes.c_int
        USE_CPP = True
        print("✓ Using C++ optimization for bonus")
except: pass

def fast_scalar_mult(k, G, curve):
    if not USE_CPP or G is None: return curve.scalar_multiply(k, G)
    rx, ry = ctypes.c_int64(), ctypes.c_int64()
    valid = ecc_lib.scalar_mult(ctypes.c_int64(k), ctypes.c_int64(G[0]), ctypes.c_int64(G[1]),
                                  ctypes.c_int64(curve.a), ctypes.c_int64(curve.b), ctypes.c_int64(curve.p),
                                  ctypes.byref(rx), ctypes.byref(ry))
    return None if valid == 0 else (rx.value, ry.value)

def fast_point_add(P, Q, curve):
    if not USE_CPP or P is None or Q is None: return curve.add(P, Q)
    rx, ry = ctypes.c_int64(), ctypes.c_int64()
    valid = ecc_lib.point_add(ctypes.c_int64(P[0]), ctypes.c_int64(P[1]),
                               ctypes.c_int64(Q[0]), ctypes.c_int64(Q[1]),
                               ctypes.c_int64(curve.a), ctypes.c_int64(curve.b), ctypes.c_int64(curve.p),
                               ctypes.c_int64(curve.n if hasattr(curve, 'n') else 0),
                               ctypes.byref(rx), ctypes.byref(ry))
    return None if valid == 0 else (rx.value, ry.value)


def pollard_rho_with_interval(curve: EllipticCurve, G: Point, Q: Point, n: int,
                              lower: int, upper: int, max_iter: int = 1000000) -> Optional[int]:
    """
    Pollard's Rho restricted to interval [lower, upper].
    Starts random walk from interval midpoint.
    """
    # Start from middle of interval
    start = (lower + upper) // 2
    
    # Initialize tortoise and hare
    a_t, b_t = start, 0
    R_t = fast_scalar_mult(a_t, G, curve)
    
    a_h, b_h = start, 0
    R_h = R_t
    
    def step(R: Point, a: int, b: int) -> tuple:
        """Single iteration step with interval awareness."""
        x = R.x if R != Point.INFINITY else 0 # type: ignore
        partition = x % 3
        
        if partition == 0:
            # Add Q
            R_new = fast_point_add(R, Q, curve)
            a_new = a
            b_new = (b + 1) % n
        elif partition == 1:
            # Double
            R_new = fast_point_add(R, R, curve)
            a_new = (2 * a) % n
            b_new = (2 * b) % n
        else:
            # Add G
            R_new = fast_point_add(R, G, curve)
            a_new = (a + 1) % n
            b_new = b
        
        # Keep a coefficient in interval if possible
        if not (lower <= a_new <= upper):
            # Adjust to stay near interval
            a_new = lower + (a_new % (upper - lower + 1))
        
        return R_new, a_new, b_new
    
    for i in range(max_iter):
        # Tortoise: one step
        R_t, a_t, b_t = step(R_t, a_t, b_t)
        
        # Hare: two steps
        R_h, a_h, b_h = step(R_h, a_h, b_h)
        R_h, a_h, b_h = step(R_h, a_h, b_h)
        
        # Check collision
        if R_t == R_h:
            r = (a_t - a_h) % n
            s = (b_h - b_t) % n
            
            if s == 0:
                # Restart with different initial point
                a_t = lower + (i % (upper - lower + 1))
                b_t = 0
                R_t = fast_scalar_mult(a_t, G, curve)
                a_h, b_h, R_h = a_t, b_t, R_t
                continue
            
            # Solve: s*d ≡ r (mod n)
            from utils import mod_inv
            try:
                s_inv = mod_inv(s, n)
                d = (r * s_inv) % n
                
                # Verify in interval
                if lower <= d <= upper:
                    test = fast_scalar_mult(d, G, curve)
                    if test == Q:
                        return d
            except:
                pass
    
    return None


def pollard_rho_with_lsb_leak(curve: EllipticCurve, G: Point, Q: Point, n: int,
                               leaked_value: int, mask: int, max_iter: int = 1000000) -> Optional[int]:
    """
    Pollard's Rho with LSB leak.
    Restricts search to candidates matching leaked bits.
    """
    step_size = mask + 1
    
    # Start from first candidate
    start = leaked_value
    
    # Try candidates matching LSB
    for offset in range(0, min(n, max_iter * step_size), step_size):
        candidate = (leaked_value + offset) % n
        
        # Quick check
        test = fast_scalar_mult(candidate, G, curve)
        if test == Q:
            return candidate
    
    return None


def main():
    """Demonstrate Pollard's Rho with partial key leakage."""
    script_dir = Path(__file__).parent
    
    if len(sys.argv) > 1:
        input_path = Path(sys.argv[1])
    else:
        input_path = script_dir / 'input' / 'test_1.txt'
    
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)
    
    # Load test case
    try:
        p, a, b, G, n, Q = load_input(input_path)
        curve = EllipticCurve(a, b, p)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    
    print("="*70)
    print("POLLARD'S RHO WITH PARTIAL KEY LEAKAGE")
    print("="*70)
    print(f"Curve: y² = x³ + {a}x + {b} (mod {p})")
    print(f"Order n = {n}")
    print(f"Task: Find d such that Q = d*G")
    print()
    
    # Note: Pollard Rho is probabilistic, so we use shorter timeouts
    max_time = 10.0  # 10 second timeout
    
    # Test 1: Bounded interval leaks
    print("="*70)
    print("SCENARIO 1: Bounded Interval Leaks")
    print("="*70)
    
    for interval_pct in [10, 5, 1]:
        interval_size = max(1000, n * interval_pct // 100)
        lower, upper = KeyLeakage.leak_bounded_interval(42, n, interval_size)
        
        # For demo, assume secret is in interval (we don't actually know it)
        # In practice, the interval would come from timing attacks, etc.
        
        print(f"\n{format_leak_info('interval', lower=lower, upper=upper, n=n)}")
        reduction, pct = calculate_search_reduction(n, upper - lower + 1)
        print(f"Search space: {upper - lower + 1:,} candidates ({pct} of original)")
        print(f"Reduction: {reduction:.2f}x")
        
        print(f"Running Pollard Rho with interval [{lower}, {upper}]...")
        start = time.perf_counter()
        
        # Try for max_time seconds
        timeout = False
        result = None
        elapsed = 0
        
        while elapsed < max_time:
            result = pollard_rho_with_interval(curve, G, Q, n, lower, upper, max_iter=10000)
            elapsed = time.perf_counter() - start
            
            if result is not None:
                break
        else:
            timeout = True
        
        if result is not None:
            print(f"✓ Found: d = {result}")
            print(f"Time: {elapsed:.6f}s")
        else:
            print(f"⊘ Timeout after {max_time}s (probabilistic method)")
    
    # Test 2: LSB leak (converts to enumeration)
    print("\n" + "="*70)
    print("SCENARIO 2: LSB Bit Leaks")
    print("="*70)
    
    for num_bits in [12, 14, 16]:
        leaked_value, mask = KeyLeakage.leak_lsb_bits(42, num_bits)
        
        print(f"\n{format_leak_info('lsb', num_bits=num_bits, leaked_value=leaked_value)}")
        candidates = candidates_from_lsb_leak(n, leaked_value, mask)
        print(f"Candidate space: {len(candidates):,} values")
        reduction, pct = calculate_search_reduction(n, len(candidates))
        print(f"Reduction: {reduction:.2f}x")
        
        print(f"Running Pollard Rho with LSB leak (uses enumeration)...")
        start = time.perf_counter()
        
        timeout = False
        result = None
        elapsed = 0
        
        # LSB leak reduces to enumeration
        while elapsed < max_time:
            result = pollard_rho_with_lsb_leak(curve, G, Q, n, leaked_value, mask, max_iter=100000)
            elapsed = time.perf_counter() - start
            
            if result is not None:
                break
        else:
            timeout = True
        
        if result is not None:
            print(f"✓ Found: d = {result}")
            print(f"Time: {elapsed:.6f}s")
        else:
            print(f"⊘ Timeout after {max_time}s")
    
    print("\n" + "="*70)
    print("KEY INSIGHT:")
    print("Pollard Rho's random walk can be restricted to leaked intervals!")
    print("However, strong leaks (LSB) make enumeration more efficient.")
    print("Interval leaks reduce collision probability space meaningfully.")
    print("="*70)


if __name__ == "__main__":
    main()
