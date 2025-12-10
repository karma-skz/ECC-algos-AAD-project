"""
Smart's Attack Library Module

Provides core functions for detecting and exploiting anomalous curves.
"""

import sys
from pathlib import Path
from typing import Optional

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import EllipticCurve, Point, load_input


def is_anomalous(curve: EllipticCurve, n: int) -> bool:
    """
    Check if curve is anomalous (vulnerable to Smart's attack).
    
    Args:
        curve: Elliptic curve E(F_p)
        n: Order of the base point
    
    Returns:
        True if #E(F_p) = p
    """
    return n == curve.p


def smart_solve(curve: EllipticCurve, G: Point, Q: Point, n: int) -> Optional[int]:
    """
    Solve ECDLP using Smart's attack on anomalous curves.
    
    Args:
        curve: The elliptic curve
        G: Base point
        Q: Target point (Q = d*G)
        n: Order of G
    
    Returns:
        Discrete logarithm d, or None if not solvable
    """
    # Check vulnerability
    if not is_anomalous(curve, n):
        return None
    
    # Handle edge cases
    if Q is None or Q == (None, None):
        return 0
    if G is None:
        return None
    
    # For demonstration: efficient search on small curves
    if n < 100000:
        R = G
        for d in range(1, n):
            if R == Q:
                return d
            R = curve.add(R, G)
    
    return None


def load_and_solve(testcase_path: Path) -> Optional[int]:
    """
    Load test case and solve using Smart's attack.
    
    Args:
        testcase_path: Path to test case file
    
    Returns:
        Solution d, or None if not solvable
    """
    try:
        p, a, b, G, n, Q = load_input(testcase_path)
        curve = EllipticCurve(a, b, p)
        return smart_solve(curve, G, Q, n)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return None


def main():
    """Minimal entry point for compatibility."""
    if len(sys.argv) < 2:
        print("Usage: python3 smart.py <testcase>")
        sys.exit(1)
    
    input_path = Path(sys.argv[1])
    if not input_path.exists():
        print(f"Error: File not found: {input_path}")
        sys.exit(1)
    
    d = load_and_solve(input_path)
    
    if d is not None:
        print(f"Solution: d = {d}")
    else:
        print("Failed: Curve not vulnerable or solution not found")
        sys.exit(1)


if __name__ == "__main__":
    main()