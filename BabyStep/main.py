"""
Baby-Step Giant-Step (BSGS) Algorithm for ECDLP

Solves the Elliptic Curve Discrete Logarithm Problem using a meet-in-the-middle approach.
Given Q = d*G, finds d by:
1. Baby steps: Compute and store j*G for j = 0, 1, ..., m-1 where m = ceil(sqrt(n))
2. Giant steps: Compute Q - i*m*G for i = 0, 1, ..., m until a match is found

Time Complexity: O(sqrt(n))
Space Complexity: O(sqrt(n))
"""

import sys
import time
import math
from pathlib import Path
from typing import Dict, Optional

# Add parent directory to path for utils import
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import EllipticCurve, Point, load_input


def bsgs_ecdlp(curve: EllipticCurve, G: Point, Q: Point, n: int) -> Optional[int]:
    """
    Solve ECDLP using Baby-Step Giant-Step algorithm.
    
    Args:
        curve: The elliptic curve
        G: Base point (generator)
        Q: Target point (Q = d*G for unknown d)
        n: Order of base point G
    
    Returns:
        The discrete logarithm d such that Q = d*G, or None if not found
    """
    if G is None:
        return None
    if Q is None:
        return 0  # 0*G = O
    if n <= 0:
        return None
    
    m = int(math.isqrt(n)) + 1  # ceil(sqrt(n))
    
    # Baby steps: Store j*G for j = 0, 1, ..., m-1
    baby_table: Dict[Point, int] = {}
    R = None  # 0*G = O
    
    for j in range(m):
        if R not in baby_table:
            baby_table[R] = j
        R = curve.add(R, G)
    
    # Compute -m*G for giant steps
    mG = curve.scalar_multiply(m, G)
    neg_mG = curve.negate(mG)
    
    # Giant steps: Check Q + i*(-m*G) for i = 0, 1, ..., m
    Gamma = Q
    for i in range(m + 1):
        if Gamma in baby_table:
            j = baby_table[Gamma]
            d = (i * m + j) % n
            return d
        Gamma = curve.add(Gamma, neg_mG)
    
    return None


def main():
    """Main entry point for BSGS ECDLP solver."""
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
    
    m = int(math.isqrt(n)) + 1
    print(f"Solving ECDLP using Baby-Step Giant-Step...")
    print(f"Curve: y^2 = x^3 + {a}x + {b} (mod {p})")
    print(f"G = ({G[0]}, {G[1]}), Q = ({Q[0]}, {Q[1]}), n = {n}")
    print(f"m = ceil(sqrt(n)) = {m}")
    
    start_time = time.perf_counter()
    d = bsgs_ecdlp(curve, G, Q, n)
    elapsed = time.perf_counter() - start_time
    
    if d is not None:
        Q_verify = curve.scalar_multiply(d, G)
        verified = (Q_verify == Q)
        
        print(f"\n{'='*50}")
        print(f"Solution: d = {d}")
        print(f"Time: {elapsed:.6f} seconds")
        print(f"Verification: {'PASSED' if verified else 'FAILED'}")
        print(f"{'='*50}")
        
        if not verified:
            sys.exit(1)
    else:
        print(f"\nNo solution found")
        print(f"Time: {elapsed:.6f} seconds")
        sys.exit(1)


if __name__ == "__main__":
    main()
