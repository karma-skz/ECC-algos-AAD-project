"""
Pohlig-Hellman Algorithm for ECDLP

Solves the Elliptic Curve Discrete Logarithm Problem when the order n
has small prime factors. Given Q = d*G where n = q1^e1 * q2^e2 * ... * qk^ek,
finds d by:
1. Solving d mod qi^ei for each prime power
2. Combining results using Chinese Remainder Theorem

Time Complexity: O(sum(ei * (log n + sqrt(qi))))
Space Complexity: O(sqrt(largest qi))
"""

import sys
import time
import math
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# Add parent directory to path for utils import
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import EllipticCurve, Point, load_input, crt_combine


def trial_factor(n: int) -> Dict[int, int]:
    """
    Factor n using trial division.
    
    Args:
        n: Number to factor
    
    Returns:
        Dictionary mapping prime factors to their exponents
    """
    factors = {}
    d = 2
    
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    
    return factors


def bsgs_small(curve: EllipticCurve, G: Point, Q: Point, order: int) -> Optional[int]:
    """
    Baby-step Giant-step for small order subgroups.
    
    Args:
        curve: The elliptic curve
        G: Base point
        Q: Target point
        order: Order of the subgroup (should be small)
    
    Returns:
        x such that x*G = Q in range [0, order-1], or None if not found
    """
    if Q is None:
        return 0
    
    m = int(math.isqrt(order)) + 1
    
    # Baby steps: Store j*G for j = 0, 1, ..., m-1
    baby_table: Dict[Point, int] = {}
    R = None
    
    for j in range(m):
        if R not in baby_table:
            baby_table[R] = j
        R = curve.add(R, G)
    
    # Giant steps
    mG = curve.scalar_multiply(m, G)
    neg_mG = curve.negate(mG)
    
    Gamma = Q
    for i in range(m + 1):
        if Gamma in baby_table:
            j = baby_table[Gamma]
            return (i * m + j) % order
        Gamma = curve.add(Gamma, neg_mG)
    
    return None


def pohlig_hellman_ecdlp(curve: EllipticCurve, G: Point, Q: Point, n: int) -> Optional[int]:
    """
    Solve ECDLP using Pohlig-Hellman algorithm.
    
    Args:
        curve: The elliptic curve
        G: Base point (generator)
        Q: Target point (Q = d*G for unknown d)
        n: Order of base point G
    
    Returns:
        The discrete logarithm d such that Q = d*G, or None if not found
    """
    # Factor n into prime powers
    factors = trial_factor(n)
    
    if not factors:
        factors = {n: 1}
    
    congruences: List[Tuple[int, int]] = []
    
    # Solve d mod q^e for each prime power
    for q, e in factors.items():
        n_i = q ** e
        h = n // n_i
        
        # Lift points to subgroup of order n_i
        G1 = curve.scalar_multiply(h, G)
        Q1 = curve.scalar_multiply(h, Q)
        
        # Use direct BSGS on the subgroup
        d_i = bsgs_small(curve, G1, Q1, n_i)
        
        if d_i is None:
            return None
        
        congruences.append((d_i, n_i))
    
    # Combine using Chinese Remainder Theorem
    d, _ = crt_combine(congruences)
    
    return d % n


def main():
    """Main entry point for Pohlig-Hellman ECDLP solver."""
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
    
    if not curve.is_on_curve(G):
        print("Error: Base point G is not on the curve", file=sys.stderr)
        sys.exit(1)
    
    if not curve.is_on_curve(Q):
        print("Error: Target point Q is not on the curve", file=sys.stderr)
        sys.exit(1)
    
    # Verify order
    nG = curve.scalar_multiply(n, G)
    if nG is not None:
        print("Warning: n*G ≠ O; provided n may not be the exact order of G")
    
    # Factor n
    factors = trial_factor(n)
    print(f"Solving ECDLP using Pohlig-Hellman...")
    print(f"Curve: y^2 = x^3 + {a}x + {b} (mod {p})")
    print(f"G = ({G[0]}, {G[1]}), Q = ({Q[0]}, {Q[1]}), n = {n}") # type: ignore
    print(f"Factorization: {factors}")
    
    start_time = time.perf_counter()
    d = pohlig_hellman_ecdlp(curve, G, Q, n)
    elapsed = time.perf_counter() - start_time
    
    if d is not None:
        Q_verify = curve.scalar_multiply(d, G)
        verified = (Q_verify == Q)
        
        # Check against answer file if it exists
        answer_path = input_path.parent / input_path.name.replace('case_', 'answer_').replace('testcase_', 'answer_')
        expected_d = None
        if answer_path.exists():
            try:
                with open(answer_path, 'r') as f:
                    expected_d = int(f.read().strip())
            except:
                pass
        
        print(f"\n{'='*50}")
        print(f"Solution: d = {d}")
        if expected_d is not None:
            print(f"Expected: d = {expected_d}")
        print(f"Time: {elapsed:.6f} seconds")
        print(f"Verification (P=d*G): {'PASSED' if verified else 'FAILED'}")
        if expected_d is not None:
            print(f"Cross-check (vs answer file): {'PASSED' if d == expected_d else 'FAILED'}")
        print(f"{'='*50}")
        
        if not verified:
            sys.exit(1)
    else:
        print(f"\nNo solution found")
        print(f"Time: {elapsed:.6f} seconds")
        sys.exit(1)


if __name__ == "__main__":
    main()
