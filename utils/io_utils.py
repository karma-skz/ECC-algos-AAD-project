"""
Input/Output utilities for reading test cases and formatting results.
"""

from pathlib import Path
from typing import Tuple
from .ecc_utils import Point


def load_input(file_path: Path) -> Tuple[int, int, int, Point, int, Point]:
    """
    Load ECDLP problem instance from file.
    
    Expected format (5 lines):
    1) p (prime modulus)
    2) a b (curve coefficients)
    3) Gx Gy (base point coordinates)
    4) n (order of base point)
    5) Qx Qy (target point coordinates)
    
    Args:
        file_path: Path to input file
    
    Returns:
        Tuple of (p, a, b, G, n, Q)
    
    Raises:
        ValueError: If file format is invalid
        FileNotFoundError: If file doesn't exist
    """
    if not file_path.exists():
        raise FileNotFoundError(f"Input file not found: {file_path}")
    
    with file_path.open('r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    if len(lines) < 5:
        raise ValueError(
            "Input file must contain 5 non-empty lines: "
            "p | a b | Gx Gy | n | Qx Qy"
        )
    
    def parse_ints(line: str):
        return list(map(int, line.split()))
    
    # Parse each line
    p = parse_ints(lines[0])[0]
    a, b = parse_ints(lines[1])
    Gx, Gy = parse_ints(lines[2])
    n = parse_ints(lines[3])[0]
    Qx, Qy = parse_ints(lines[4])
    
    G = (Gx, Gy)
    Q = (Qx, Qy)
    
    return p, a, b, G, n, Q


def format_output(algorithm: str, p: int, a: int, b: int, G: Point, Q: Point, 
                 n: int, d: int, elapsed: float, verified: bool) -> str:
    """
    Format the output of an ECDLP solution.
    
    Args:
        algorithm: Name of the algorithm used
        p, a, b: Curve parameters
        G: Base point
        Q: Target point
        n: Order of base point
        d: Solution (discrete log)
        elapsed: Time elapsed in seconds
        verified: Whether solution was verified
    
    Returns:
        Formatted output string
    """
    if G is None or Q is None:
        raise ValueError("G and Q cannot be point at infinity")
    
    output = []
    output.append(f"Algorithm: {algorithm}")
    output.append(f"Curve: y^2 = x^3 + {a}x + {b} (mod {p})")
    output.append(f"Base point G: ({G[0]}, {G[1]})")
    output.append(f"Target point Q: ({Q[0]}, {Q[1]})")
    output.append(f"Order n: {n}")
    output.append(f"Solution d: {d}")
    output.append(f"Time: {elapsed:.6f} seconds")
    output.append(f"Verification: {'PASSED' if verified else 'FAILED'}")
    
    return '\n'.join(output)
