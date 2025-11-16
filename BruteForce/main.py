"""
Brute Force Algorithm for ECDLP

Solves the Elliptic Curve Discrete Logarithm Problem by exhaustive search.
Given Q = d*G, finds d by computing k*G for k = 1, 2, ..., n-1 until Q is found.

Time Complexity: O(n) where n is the order of the base point
Space Complexity: O(1)
"""

import sys
import time
from pathlib import Path

# Add parent directory to path for utils import
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import EllipticCurve, Point, load_input


def brute_force_ecdlp(curve: EllipticCurve, G: Point, Q: Point, n: int) -> int:
    """
    Solve ECDLP using brute force exhaustive search.
    
    Args:
        curve: The elliptic curve
        G: Base point (generator)
        Q: Target point (Q = d*G for unknown d)
        n: Order of base point G
    
    Returns:
        The discrete logarithm d such that Q = d*G
    
    Raises:
        ValueError: If no solution found
    """
    if Q is None:
        raise ValueError("Q cannot be the point at infinity")
    
    # Check k*G for k = 1, 2, ..., n-1
    R = G
    for k in range(1, n):
        if R == Q:
            return k
        R = curve.add(R, G)
    
    # Check k = n (should give point at infinity)
    if R == Q:
        return n
    
    raise ValueError("No solution found in range [1, n]")


def main():
    """Main entry point for brute force ECDLP solver."""
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
    
    print(f"Solving ECDLP using Brute Force...")
    print(f"Curve: y^2 = x^3 + {a}x + {b} (mod {p})")
    print(f"G = ({G[0]}, {G[1]}), Q = ({Q[0]}, {Q[1]}), n = {n}") # type: ignore
    
    start_time = time.perf_counter()
    try:
        d = brute_force_ecdlp(curve, G, Q, n)
        elapsed = time.perf_counter() - start_time
        
        Q_verify = curve.scalar_multiply(d, G)
        verified = (Q_verify == Q)
        
        print(f"\n{'='*50}")
        print(f"Solution: d = {d}")
        print(f"Time: {elapsed:.6f} seconds")
        print(f"Verification: {'PASSED' if verified else 'FAILED'}")
        print(f"{'='*50}")
        
        if not verified:
            sys.exit(1)
            
    except ValueError as e:
        elapsed = time.perf_counter() - start_time
        print(f"\nNo solution found: {e}")
        print(f"Time: {elapsed:.6f} seconds")
        sys.exit(1)


if __name__ == "__main__":
    main()
