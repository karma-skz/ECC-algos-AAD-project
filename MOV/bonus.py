"""
MOV Attack with Partial Key Leakage - BONUS Implementation

Demonstrates how bit leakage can reduce the DLP search space in the extension field.
Even though MOV maps ECDLP to field DLP, leaked bits still constrain the search.

Scenarios:
1. Known MSB bits: Reduces field DLP search range
2. Known LSB bits: Allows verification of candidates
3. Bounded interval: Restricts field element search
"""

import sys
import time
import ctypes
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import (EllipticCurve, Point, load_input, KeyLeakage,
                   candidates_from_msb_leak, candidates_from_lsb_leak,
                   format_leak_info, calculate_speedup, calculate_search_reduction)
from MOV.main import (find_embedding_degree, generate_irreducible_poly,
                      Polynomial, FieldElement, bsgs_field)

# Load C++ library
USE_CPP = False
try:
    lib_path = Path(__file__).parent / 'ecc_fast.so'
    if lib_path.exists():
        ecc_lib = ctypes.CDLL(str(lib_path))
        ecc_lib.scalar_mult.argtypes = [ctypes.c_int64]*6 + [ctypes.POINTER(ctypes.c_int64)]*2
        ecc_lib.scalar_mult.restype = ctypes.c_int
        USE_CPP = True
        print("✓ Using C++ optimization for bonus")
except: pass


def fast_scalar_mult(k: int, G: Point, curve: EllipticCurve) -> Point:
    """Fast scalar multiplication using C++ if available."""
    if not USE_CPP or G is None:
        return curve.scalar_multiply(k, G)
    
    gx, gy = G
    result_x = ctypes.c_int64()
    result_y = ctypes.c_int64()
    
    valid = ecc_lib.scalar_mult( # type: ignore
        ctypes.c_int64(k),
        ctypes.c_int64(gx),
        ctypes.c_int64(gy),
        ctypes.c_int64(curve.a),
        ctypes.c_int64(curve.b),
        ctypes.c_int64(curve.p),
        ctypes.byref(result_x),
        ctypes.byref(result_y)
    )
    
    if valid == 0:
        return None
    return (result_x.value, result_y.value)


def bsgs_field_with_leak(alpha: FieldElement, beta: FieldElement, n: int,
                         lower: int, upper: int) -> Optional[int]:
    """
    BSGS in extension field restricted to interval [lower, upper].
    This is the field DLP solver used after pairing computation.
    """
    from math import isqrt
    
    # Effective range
    range_size = upper - lower + 1
    m = isqrt(range_size) + 1
    
    # Adjust for lower bound: beta' = beta * alpha^(-lower)
    q = alpha.p ** alpha.mod_poly.degree()
    alpha_inv = alpha ** (q - 2)  # Fermat's little theorem
    lower_factor = alpha_inv ** lower
    beta_adjusted = beta * lower_factor
    
    # Baby steps
    baby_table = {}
    gamma = FieldElement.from_int(1, alpha.p, alpha.mod_poly)
    
    for j in range(m):
        key = tuple(gamma.poly.coeffs)
        if key not in baby_table:
            baby_table[key] = j
        gamma = gamma * alpha
    
    # Giant steps
    factor = alpha_inv ** m
    gamma = beta_adjusted
    
    for i in range(m + 1):
        key = tuple(gamma.poly.coeffs)
        if key in baby_table:
            j = baby_table[key]
            d_offset = i * m + j
            if d_offset <= range_size:
                return lower + d_offset
        gamma = gamma * factor
    
    return None


def mov_with_msb_leak(curve: EllipticCurve, G: Point, Q: Point, n: int,
                      leaked_msb: int, shift: int, max_k: int = 12) -> Optional[int]:
    """
    MOV attack with known MSB bits.
    The pairing computation is the same, but field DLP search is restricted.
    """
    p = curve.p
    
    # Find embedding degree
    k = find_embedding_degree(p, n, max_k)
    if k is None:
        return None
    
    # Generate extension field
    mod_poly = generate_irreducible_poly(k, p)
    
    # In a full implementation, we would:
    # 1. Compute pairings α = e(G, R), β = e(Q, R)
    # 2. Use leaked bits to restrict BSGS search in extension field
    
    # Since we don't have full pairing, demonstrate the concept
    print(f"  Embedding degree k = {k}")
    print(f"  Extension field: GF({p}^{k})")
    
    # Demonstrate reduced search space
    lower, upper = candidates_from_msb_leak(leaked_msb, shift, n)
    print(f"  Search range: [{lower}, {upper}] (size: {upper - lower + 1})")
    
    # In practice: d = bsgs_field_with_leak(alpha, beta, n, lower, upper)
    
    return None


def mov_with_interval_leak(curve: EllipticCurve, G: Point, Q: Point, n: int,
                            lower: int, upper: int, max_k: int = 12) -> Optional[int]:
    """MOV attack with bounded interval leak."""
    p = curve.p
    
    k = find_embedding_degree(p, n, max_k)
    if k is None:
        return None
    
    mod_poly = generate_irreducible_poly(k, p)
    
    print(f"  Embedding degree k = {k}")
    print(f"  Extension field: GF({p}^{k})")
    print(f"  Interval: [{lower}, {upper}] (size: {upper - lower + 1})")
    
    # In practice: d = bsgs_field_with_leak(alpha, beta, n, lower, upper)
    
    return None


def main():
    """Demonstrate MOV attack with partial key leakage."""
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
    
    # Load actual answer if available
    answer_file = input_path.parent / f"answer_{input_path.stem.split('_')[-1]}.txt"
    d_actual = None
    if answer_file.exists():
        with answer_file.open('r') as f:
            d_actual = int(f.read().strip())
    
    total_bits = n.bit_length()
    
    print("="*70)
    print("MOV ATTACK WITH PARTIAL KEY LEAKAGE")
    print("="*70)
    print(f"Curve: y² = x³ + {a}x + {b} (mod {p})")
    print(f"Order n = {n} ({total_bits} bits)")
    if d_actual:
        print(f"Actual secret: d = {d_actual}")
    print()
    
    # Scenario 1: Known MSB bits
    print("="*70)
    print("SCENARIO 1: Known MSB Bits")
    print("="*70)
    print()
    
    if d_actual:
        for num_bits in [8, 12, 16]:
            leaked_msb, shift = KeyLeakage.leak_msb_bits(d_actual, num_bits, total_bits)
            lower, upper = candidates_from_msb_leak(leaked_msb, shift, n)
            search_space = upper - lower + 1
            
            print(f"{format_leak_info('msb', num_bits=num_bits, leaked=leaked_msb)}")
            print(f"Search space: {search_space:,} candidates (vs {n:,} original)")
            reduction, pct = calculate_search_reduction(n, search_space)
            print(f"Reduction: {reduction:.2f}x ({pct} smaller)")
            
            start = time.perf_counter()
            d = mov_with_msb_leak(curve, G, Q, n, leaked_msb, shift, max_k=8)
            elapsed = time.perf_counter() - start
            
            print(f"Time: {elapsed:.6f}s")
            print(f"⚠ Full pairing not implemented")
            print()
    else:
        print("⚠ No answer file available for demonstration")
        print()
    
    # Scenario 2: Bounded interval
    print("="*70)
    print("SCENARIO 2: Bounded Interval")
    print("="*70)
    print()
    
    if d_actual:
        for leak_pct in [0.10, 0.05, 0.01]:
            lower, upper = KeyLeakage.leak_bounded_interval(d_actual, n, leak_pct)
            search_space = upper - lower + 1
            
            print(f"{format_leak_info('interval', lower=lower, upper=upper)}")
            print(f"Leak percentage: {leak_pct*100:.1f}% of n")
            reduction, pct = calculate_search_reduction(n, search_space)
            print(f"Reduction: {reduction:.2f}x ({pct} smaller)")
            
            start = time.perf_counter()
            d = mov_with_interval_leak(curve, G, Q, n, lower, upper, max_k=8)
            elapsed = time.perf_counter() - start
            
            print(f"Time: {elapsed:.6f}s")
            print(f"⚠ Full pairing not implemented")
            print()
    else:
        print("⚠ No answer file available for demonstration")
        print()
    
    # Educational summary
    print("="*70)
    print("KEY INSIGHTS: MOV Attack with Leakage")
    print("="*70)
    print()
    print("1. Leakage Impact:")
    print("   • MSB leaks reduce the DLP search space in F_p^k")
    print("   • Works because pairing preserves linear relationships")
    print("   • If Q = d*G, then e(Q,R) = e(G,R)^d")
    print()
    print("2. Advantages over direct ECDLP:")
    print("   • Can use field-specific optimizations")
    print("   • Index calculus possible in some extension fields")
    print("   • Parallel search easier in field representation")
    print()
    print("3. Practical Considerations:")
    print("   • Only works when embedding degree k is small")
    print("   • Modern secure curves have k > 10^6 (infeasible)")
    print("   • Pairing computation itself is expensive")
    print("   • Full implementation requires Miller's algorithm")
    print()
    print("4. With Leakage:")
    print("   • MSB leak: Reduces field BSGS search proportionally")
    print("   • LSB leak: Allows candidate verification")
    print("   • Interval leak: Restricts search to subset of F_p^k*")
    print("="*70)


if __name__ == "__main__":
    main()
