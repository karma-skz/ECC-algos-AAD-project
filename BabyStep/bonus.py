"""
Baby-Step Giant-Step with Partial Key Leakage - BONUS Implementation (C++ Optimized)

Demonstrates how known bits or reduced range optimize BSGS.
Uses C++ backend for fast operations.
"""

import sys
import time
import math
import ctypes
from pathlib import Path
from typing import Optional, Dict

sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import (EllipticCurve, Point, load_input, KeyLeakage,
                   candidates_from_msb_leak, format_leak_info,
                   calculate_speedup, calculate_search_reduction)

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


def bsgs_with_interval(curve: EllipticCurve, G: Point, Q: Point,
                      lower: int, upper: int) -> Optional[int]:
    """
    BSGS within a bounded interval [lower, upper].
    Uses effective range size instead of full order n.
    """
    if G is None or Q is None:
        return None
    
    # Effective range
    range_size = upper - lower + 1
    m = int(math.isqrt(range_size)) + 1
    
    # Adjust Q to account for lower bound: Q' = Q - lower*G
    G_lower = fast_scalar_mult(lower, G, curve)
    Q_adjusted = fast_point_add(Q, curve.negate(G_lower), curve)
    
    # Baby steps: j*G for j = 0, 1, ..., m-1
    baby_table: Dict[Point, int] = {}
    R = None
    for j in range(m):
        if R not in baby_table:
            baby_table[R] = j
        R = fast_point_add(R, G, curve)
    
    # Giant steps: Q' - i*m*G
    mG = fast_scalar_mult(m, G, curve)
    neg_mG = curve.negate(mG)
    
    Gamma = Q_adjusted
    for i in range(m + 1):
        if Gamma in baby_table:
            j = baby_table[Gamma]
            d_offset = i * m + j
            if d_offset <= range_size:
                return lower + d_offset
        Gamma = fast_point_add(Gamma, neg_mG, curve)
    
    return None


def bsgs_with_msb_leak(curve: EllipticCurve, G: Point, Q: Point, n: int,
                      leaked_msb: int, shift: int) -> Optional[int]:
    """
    BSGS with known MSB bits.
    Reduces search to range defined by MSB.
    """
    lower, upper = candidates_from_msb_leak(leaked_msb, shift, n)
    return bsgs_with_interval(curve, G, Q, lower, upper)


def main():
    """Demonstrate BSGS with partial key leakage."""
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
    
    # Load actual answer
    answer_file = input_path.parent / f"answer_{input_path.stem.split('_')[1]}.txt"
    if not answer_file.exists():
        print("Error: Answer file not found", file=sys.stderr)
        sys.exit(1)
    
    with answer_file.open('r') as f:
        d_actual = int(f.read().strip())
    
    total_bits = n.bit_length()
    
    print("="*70)
    print("BABY-STEP GIANT-STEP WITH PARTIAL KEY LEAKAGE")
    print("="*70)
    print(f"Curve: y² = x³ + {a}x + {b} (mod {p})")
    print(f"Order n = {n} ({total_bits} bits)")
    print(f"Actual secret: d = {d_actual}")
    print()
    
    # Test 1: Known MSB bits
    print("="*70)
    print("SCENARIO 1: Known MSB Bits")
    print("="*70)
    
    for num_bits in [8, 12, 16]:
        leaked_msb, shift = KeyLeakage.leak_msb_bits(d_actual, num_bits, total_bits)
        lower, upper = candidates_from_msb_leak(leaked_msb, shift, n)
        search_space = upper - lower + 1
        
        print(f"\n{format_leak_info('msb', num_bits=num_bits, leaked=leaked_msb)}")
        print(f"Search range: [{lower:,}, {upper:,}]")
        print(f"Search space: {search_space:,} candidates (vs {n:,} original)")
        reduction, pct = calculate_search_reduction(n, search_space)
        print(f"Reduction: {reduction:.2f}x ({pct} smaller)")
        print(f"BSGS m = {int(math.isqrt(search_space)) + 1:,} (vs {int(math.isqrt(n)) + 1:,} original)")
        
        start = time.perf_counter()
        d = bsgs_with_msb_leak(curve, G, Q, n, leaked_msb, shift)
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
        print(f"BSGS m = {int(math.isqrt(search_space)) + 1:,} (vs {int(math.isqrt(n)) + 1:,} original)")
        
        start = time.perf_counter()
        d = bsgs_with_interval(curve, G, Q, lower, upper)
        elapsed = time.perf_counter() - start
        
        if d == d_actual:
            print(f"✓ Found: d = {d}")
            print(f"Time: {elapsed:.6f}s")
        else:
            print(f"✗ Failed or incorrect result")
    
    # Comparison with standard BSGS
    print("\n" + "="*70)
    print("COMPARISON: Standard vs Leaked")
    print("="*70)
    
    print("\nRunning standard BSGS (no leak)...")
    from BabyStep.main import bsgs_ecdlp
    start = time.perf_counter()
    d_standard = bsgs_ecdlp(curve, G, Q, n)
    time_standard = time.perf_counter() - start
    print(f"Time: {time_standard:.6f}s")
    
    # Best leak scenario
    leaked_msb, shift = KeyLeakage.leak_msb_bits(d_actual, 16, total_bits)
    print(f"\nRunning with 16-bit MSB leak...")
    start = time.perf_counter()
    d_leak = bsgs_with_msb_leak(curve, G, Q, n, leaked_msb, shift)
    time_leak = time.perf_counter() - start
    print(f"Time: {time_leak:.6f}s")
    
    speedup = calculate_speedup(time_standard, time_leak)
    print(f"\n→ Speedup: {speedup:.2f}x faster with leak!")
    print("="*70)


if __name__ == "__main__":
    main()
