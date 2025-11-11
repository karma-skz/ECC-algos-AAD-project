"""
Brute Force with Partial Key Leakage - BONUS Implementation (C++ Optimized)

Demonstrates how known bits or bounded intervals drastically reduce search space.
Uses C++ backend for fast scalar multiplication.
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
        ecc_lib.point_add.argtypes = [ctypes.c_int64]*6 + [ctypes.POINTER(ctypes.c_int64)]*2
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
    if not USE_CPP or P is None: return curve.add(P, Q)
    if Q is None: return P
    rx, ry = ctypes.c_int64(), ctypes.c_int64()
    valid = ecc_lib.point_add(ctypes.c_int64(P[0]), ctypes.c_int64(P[1]),
                               ctypes.c_int64(Q[0]), ctypes.c_int64(Q[1]),
                               ctypes.c_int64(curve.a), ctypes.c_int64(curve.p),
                               ctypes.byref(rx), ctypes.byref(ry))
    return None if valid == 0 else (rx.value, ry.value)


def brute_force_with_lsb_leak(curve: EllipticCurve, G: Point, Q: Point, n: int,
                               leaked_lsb: int, mask: int) -> Optional[int]:
    """
    Brute force search with known LSB bits.
    Only checks candidates consistent with leaked LSB.
    """
    step = mask + 1
    R = fast_scalar_mult(leaked_lsb, G, curve)
    step_G = fast_scalar_mult(step, G, curve)
    
    k = leaked_lsb
    while k < n:
        if R == Q:
            return k
        R = fast_point_add(R, step_G, curve)
        k += step
    
    return None


def brute_force_with_interval(curve: EllipticCurve, G: Point, Q: Point, n: int,
                              lower: int, upper: int) -> Optional[int]:
    """
    Brute force search within bounded interval.
    """
    R = fast_scalar_mult(lower, G, curve)
    
    for k in range(lower, upper + 1):
        if R == Q:
            return k
        R = fast_point_add(R, G, curve)
    
    return None


def main():
    """Demonstrate brute force with partial key leakage."""
    script_dir = Path(__file__).parent
    
    if len(sys.argv) > 1:
        input_path = Path(sys.argv[1])
    else:
        input_path = script_dir.parent / 'input' / 'testcase_1.txt'
    
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
    
    # Load actual answer for leak generation
    answer_file = input_path.parent / f"answer_{input_path.stem.split('_')[1]}.txt"
    if not answer_file.exists():
        print("Error: Answer file not found", file=sys.stderr)
        sys.exit(1)
    
    with answer_file.open('r') as f:
        d_actual = int(f.read().strip())
    
    print("="*70)
    print("BRUTE FORCE WITH PARTIAL KEY LEAKAGE")
    print("="*70)
    print(f"Curve: y² = x³ + {a}x + {b} (mod {p})")
    print(f"Order n = {n}")
    print(f"Actual secret: d = {d_actual}")
    print()
    
    # Test 1: Known LSB bits
    print("="*70)
    print("SCENARIO 1: Known LSB Bits")
    print("="*70)
    
    for num_bits in [8, 12, 16]:
        leaked_lsb, mask = KeyLeakage.leak_lsb_bits(d_actual, num_bits)
        search_space = (n + mask) // (mask + 1)
        
        print(f"\n{format_leak_info('lsb', num_bits=num_bits, leaked=leaked_lsb)}")
        print(f"Search space: {search_space:,} candidates (vs {n:,} original)")
        reduction, pct = calculate_search_reduction(n, search_space)
        print(f"Reduction: {reduction:.2f}x ({pct} smaller)")
        
        start = time.perf_counter()
        d = brute_force_with_lsb_leak(curve, G, Q, n, leaked_lsb, mask)
        elapsed = time.perf_counter() - start
        
        if d == d_actual:
            print(f"✓ Found: d = {d}")
            print(f"Time: {elapsed:.6f}s")
        else:
            print(f"✗ Failed or incorrect result")
    
    # Test 2: Bounded interval
    print("\n" + "="*70)
    print("SCENARIO 2: Bounded Interval")
    print("="*70)
    
    for leak_pct in [0.10, 0.05, 0.01]:
        lower, upper = KeyLeakage.leak_bounded_interval(d_actual, n, leak_pct)
        search_space = upper - lower + 1
        
        print(f"\n{format_leak_info('interval', lower=lower, upper=upper)}")
        print(f"Leak percentage: {leak_pct*100:.1f}% of n")
        reduction, pct = calculate_search_reduction(n, search_space)
        print(f"Reduction: {reduction:.2f}x ({pct} smaller)")
        
        start = time.perf_counter()
        d = brute_force_with_interval(curve, G, Q, n, lower, upper)
        elapsed = time.perf_counter() - start
        
        if d == d_actual:
            print(f"✓ Found: d = {d}")
            print(f"Time: {elapsed:.6f}s")
        else:
            print(f"✗ Failed or incorrect result")
    
    # Comparison with standard brute force
    print("\n" + "="*70)
    print("COMPARISON: Standard vs Leaked")
    print("="*70)
    
    print("\nRunning standard brute force (no leak)...")
    from BruteForce.main import brute_force_ecdlp
    start = time.perf_counter()
    d_standard = brute_force_ecdlp(curve, G, Q, n)
    time_standard = time.perf_counter() - start
    print(f"Time: {time_standard:.6f}s")
    
    # Best leak scenario
    leaked_lsb, mask = KeyLeakage.leak_lsb_bits(d_actual, 16)
    print(f"\nRunning with 16-bit LSB leak...")
    start = time.perf_counter()
    d_leak = brute_force_with_lsb_leak(curve, G, Q, n, leaked_lsb, mask)
    time_leak = time.perf_counter() - start
    print(f"Time: {time_leak:.6f}s")
    
    speedup = calculate_speedup(time_standard, time_leak)
    print(f"\n→ Speedup: {speedup:.2f}x faster with leak!")
    print("="*70)


if __name__ == "__main__":
    main()
