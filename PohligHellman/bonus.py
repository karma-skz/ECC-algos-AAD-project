"""
Pohlig-Hellman with Partial Key Leakage - BONUS Implementation (C++ Optimized)

Demonstrates how residue leaks eliminate factors from CRT reconstruction.
Uses C++ backend for fast operations.
"""

import sys
import time
import ctypes
from pathlib import Path
from typing import Optional, List, Tuple, Dict

sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import (EllipticCurve, Point, load_input, KeyLeakage, crt_combine,
                   format_leak_info, calculate_speedup, calculate_search_reduction)
from PohligHellman.main import trial_factor, bsgs_small

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


def pohlig_hellman_with_residue_leak(curve: EllipticCurve, G: Point, Q: Point, n: int,
                                    leaked_residues: List[Tuple[int, int]]) -> Optional[int]:
    """
    Pohlig-Hellman with known residues.
    Skips computation for already-known factors.
    """
    # Factor n
    factors = trial_factor(n)
    
    if not factors:
        factors = {n: 1}
    
    congruences: List[Tuple[int, int]] = []
    
    # Add leaked residues directly
    leaked_moduli = {m for _, m in leaked_residues}
    for residue, modulus in leaked_residues:
        congruences.append((residue, modulus))
    
    # Solve only for non-leaked factors
    for q, e in factors.items():
        n_i = q ** e
        
        # Skip if we already have this residue
        if n_i in leaked_moduli:
            continue
        
        h = n // n_i
        
        # Lift points to subgroup
        G1 = fast_scalar_mult(h, G, curve)
        Q1 = fast_scalar_mult(h, Q, curve)
        
        # Solve in subgroup
        d_i = bsgs_small(curve, G1, Q1, n_i)
        
        if d_i is None:
            return None
        
        congruences.append((d_i, n_i))
    
    # Combine using CRT
    d, _ = crt_combine(congruences)
    
    return d % n


def pohlig_hellman_with_partial_residues(curve: EllipticCurve, G: Point, Q: Point, n: int,
                                         num_leaked: int) -> Tuple[Optional[int], List[Tuple[int, int]]]:
    """
    Simulate partial residue leak and solve.
    
    Returns:
        (solution, leaked_residues)
    """
    # Factor n to find small primes
    factors = trial_factor(n)
    if not factors:
        return None, []
    
    # Select smallest primes to leak
    sorted_factors = sorted(factors.items(), key=lambda x: x[0] ** x[1])
    leaked_factors = sorted_factors[:num_leaked]
    
    # Generate leaked residues (in practice, these would be obtained through side-channels)
    # Here we compute them from the actual solution for demonstration
    h = n // n  # We don't actually know d yet, so this is just for structure
    
    # For demo, we'll leak the factors but still solve the problem
    # In practice, you'd get the residues from side-channel attacks
    
    return None, []


def main():
    """Demonstrate Pohlig-Hellman with partial key leakage."""
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
    
    # Factor n to understand structure
    factors = trial_factor(n)
    
    print("="*70)
    print("POHLIG-HELLMAN WITH PARTIAL KEY LEAKAGE")
    print("="*70)
    print(f"Curve: y² = x³ + {a}x + {b} (mod {p})")
    print(f"Order n = {n}")
    print(f"Factorization: {factors}")
    print(f"Actual secret: d = {d_actual}")
    print()
    
    # Test: Known residues for small primes
    print("="*70)
    print("SCENARIO: Known Residues for Small Primes")
    print("="*70)
    
    # Sort factors by size
    sorted_factors = sorted(factors.items(), key=lambda x: x[0] ** x[1])
    
    for num_leaked in [1, 2, 3]:
        if num_leaked > len(factors):
            break
        
        # Leak residues for smallest primes
        leaked_moduli = [q ** e for q, e in sorted_factors[:num_leaked]]
        leaked_residues = KeyLeakage.leak_residues(d_actual, leaked_moduli)
        
        # Calculate work saved
        total_work = sum(q ** e for q, e in factors.items())
        leaked_work = sum(leaked_moduli)
        remaining_work = total_work - leaked_work
        
        print(f"\n{format_leak_info('residues', residues=leaked_residues)}")
        print(f"Factors leaked: {num_leaked}/{len(factors)}")
        print(f"Work saved: {leaked_work:,} (remaining: {remaining_work:,})")
        reduction, pct = calculate_search_reduction(total_work, remaining_work)
        print(f"Complexity reduction: {reduction:.2f}x ({pct} easier)")
        
        start = time.perf_counter()
        d = pohlig_hellman_with_residue_leak(curve, G, Q, n, leaked_residues)
        elapsed = time.perf_counter() - start
        
        if d == d_actual:
            print(f"✓ Found: d = {d}")
            print(f"Time: {elapsed:.6f}s")
        else:
            print(f"✗ Failed or incorrect result: {d}")
    
    # Comparison with standard Pohlig-Hellman
    print("\n" + "="*70)
    print("COMPARISON: Standard vs Leaked")
    print("="*70)
    
    print("\nRunning standard Pohlig-Hellman (no leak)...")
    from PohligHellman.main import pohlig_hellman_ecdlp
    start = time.perf_counter()
    d_standard = pohlig_hellman_ecdlp(curve, G, Q, n)
    time_standard = time.perf_counter() - start
    print(f"Time: {time_standard:.6f}s")
    
    # Best leak scenario
    if len(sorted_factors) >= 2:
        leaked_moduli = [q ** e for q, e in sorted_factors[:2]]
        leaked_residues = KeyLeakage.leak_residues(d_actual, leaked_moduli)
        
        print(f"\nRunning with {len(leaked_residues)} residue leaks...")
        start = time.perf_counter()
        d_leak = pohlig_hellman_with_residue_leak(curve, G, Q, n, leaked_residues)
        time_leak = time.perf_counter() - start
        print(f"Time: {time_leak:.6f}s")
        
        speedup = calculate_speedup(time_standard, time_leak)
        print(f"\n→ Speedup: {speedup:.2f}x faster with leak!")
    
    print("="*70)
    print("\nKEY INSIGHT:")
    print("Pohlig-Hellman's strength becomes its weakness with residue leaks!")
    print("Each leaked residue eliminates entire subproblems from the CRT system.")
    print("="*70)


if __name__ == "__main__":
    main()
