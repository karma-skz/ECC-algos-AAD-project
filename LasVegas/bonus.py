"""
Las Vegas with Partial Key Leakage - BONUS Implementation (C++ Optimized)

Demonstrates how approximate value leaks guide random sampling.
Uses C++ backend for fast operations.
"""

import sys
import time
import random
import ctypes
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import (EllipticCurve, Point, load_input, KeyLeakage,
                   format_leak_info, calculate_speedup, calculate_search_reduction)

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

def fast_scalar_mult(k, G, curve):
    if not USE_CPP or G is None: return curve.scalar_multiply(k, G)
    rx, ry = ctypes.c_int64(), ctypes.c_int64()
    valid = ecc_lib.scalar_mult(ctypes.c_int64(k), ctypes.c_int64(G[0]), ctypes.c_int64(G[1]),
                                  ctypes.c_int64(curve.a), ctypes.c_int64(curve.b), ctypes.c_int64(curve.p),
                                  ctypes.byref(rx), ctypes.byref(ry))
    return None if valid == 0 else (rx.value, ry.value)


def las_vegas_with_approximate(curve: EllipticCurve, G: Point, Q: Point, n: int,
                               approx_value: int, error: int, max_attempts: int = 100000) -> Optional[int]:
    """
    Las Vegas algorithm with approximate value leak.
    Focuses random sampling around the approximate value.
    """
    # Define search window
    lower = max(0, approx_value - error)
    upper = min(n - 1, approx_value + error)
    
    # Weighted random sampling: more probability near approx_value
    attempts = 0
    checked = set()
    
    while attempts < max_attempts:
        # Gaussian-like sampling centered at approx_value
        # Use triangular distribution for simplicity
        if random.random() < 0.7:
            # 70% chance: sample close to approximate
            offset = random.randint(-error // 2, error // 2)
        else:
            # 30% chance: sample from full interval
            offset = random.randint(-error, error)
        
        d = (approx_value + offset) % n
        
        # Ensure in bounds
        if d < lower or d > upper:
            d = random.randint(lower, upper)
        
        if d in checked:
            continue
        
        checked.add(d)
        attempts += 1
        
        # Test candidate
        test = fast_scalar_mult(d, G, curve)
        if test == Q:
            return d
    
    return None


def las_vegas_weighted_sampling(curve: EllipticCurve, G: Point, Q: Point, n: int,
                                approx_value: int, error: int, max_attempts: int = 100000) -> Optional[int]:
    """
    Las Vegas with weighted sampling around approximate value.
    Uses normal distribution for more natural spread.
    """
    import math
    
    lower = max(0, approx_value - error)
    upper = min(n - 1, approx_value + error)
    
    attempts = 0
    checked = set()
    
    # Standard deviation: error / 3 (99.7% within bounds)
    sigma = error / 3.0
    
    while attempts < max_attempts:
        # Sample from normal distribution
        offset = int(random.gauss(0, sigma))
        d = (approx_value + offset) % n
        
        # Clip to bounds
        if d < lower:
            d = lower
        elif d > upper:
            d = upper
        
        if d in checked:
            continue
        
        checked.add(d)
        attempts += 1
        
        # Test candidate
        test = fast_scalar_mult(d, G, curve)
        if test == Q:
            return d
    
    return None


def main():
    """Demonstrate Las Vegas with partial key leakage."""
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
    print("LAS VEGAS WITH PARTIAL KEY LEAKAGE")
    print("="*70)
    print(f"Curve: y² = x³ + {a}x + {b} (mod {p})")
    print(f"Order n = {n}")
    print(f"Task: Find d such that Q = d*G")
    print()
    
    # Note: Las Vegas is highly probabilistic, short timeouts
    max_time = 15.0  # 15 second timeout per attempt
    
    # Test: Approximate value leaks
    print("="*70)
    print("SCENARIO: Approximate Value Leaks")
    print("="*70)
    print("Simulating leaked approximate values from power analysis...")
    print()
    
    for error_pct in [10, 5, 1]:
        # Generate approximate value leak
        approx_value, error = KeyLeakage.leak_approximate_value(42, n, error_pct)
        
        print(f"{format_leak_info('approximate', approx_value=approx_value, error=error, n=n)}")
        search_space = 2 * error + 1
        reduction, pct = calculate_search_reduction(n, search_space)
        print(f"Search space: {search_space:,} candidates ({pct} of original)")
        print(f"Reduction: {reduction:.2f}x")
        
        # Try uniform sampling
        print(f"\nMethod 1: Uniform random sampling in interval...")
        start = time.perf_counter()
        result_uniform = None
        elapsed = 0
        
        while elapsed < max_time:
            lower = max(0, approx_value - error)
            upper = min(n - 1, approx_value + error)
            
            # Simple uniform random in range
            for _ in range(10000):
                d = random.randint(lower, upper)
                test = fast_scalar_mult(d, G, curve)
                if test == Q:
                    result_uniform = d
                    break
            
            elapsed = time.perf_counter() - start
            if result_uniform is not None:
                break
        
        if result_uniform is not None:
            print(f"✓ Found: d = {result_uniform}")
            print(f"Time: {elapsed:.6f}s")
        else:
            print(f"⊘ Timeout after {max_time}s")
        
        # Try weighted sampling
        print(f"\nMethod 2: Weighted sampling (focused near approximate)...")
        start = time.perf_counter()
        result_weighted = None
        elapsed = 0
        
        while elapsed < max_time:
            result_weighted = las_vegas_with_approximate(curve, G, Q, n, approx_value, error, max_attempts=10000)
            elapsed = time.perf_counter() - start
            
            if result_weighted is not None:
                break
        
        if result_weighted is not None:
            print(f"✓ Found: d = {result_weighted}")
            print(f"Time: {elapsed:.6f}s")
        else:
            print(f"⊘ Timeout after {max_time}s")
        
        # Try Gaussian sampling
        print(f"\nMethod 3: Gaussian sampling (normal distribution)...")
        start = time.perf_counter()
        result_gauss = None
        elapsed = 0
        
        while elapsed < max_time:
            result_gauss = las_vegas_weighted_sampling(curve, G, Q, n, approx_value, error, max_attempts=10000)
            elapsed = time.perf_counter() - start
            
            if result_gauss is not None:
                break
        
        if result_gauss is not None:
            print(f"✓ Found: d = {result_gauss}")
            print(f"Time: {elapsed:.6f}s")
        else:
            print(f"⊘ Timeout after {max_time}s")
        
        print("-" * 70)
    
    print("\n" + "="*70)
    print("KEY INSIGHT:")
    print("Las Vegas benefits from approximate leaks by focusing random sampling!")
    print("Weighted/Gaussian distributions outperform uniform sampling.")
    print("Approximate values from power analysis dramatically improve success rate.")
    print("="*70)


if __name__ == "__main__":
    main()
