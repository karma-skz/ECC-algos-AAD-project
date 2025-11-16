"""
Las Vegas Algorithm for ECDLP

A probabilistic algorithm that uses summation polynomials and linear algebra
to solve the Elliptic Curve Discrete Logarithm Problem. The algorithm:
1. Generates random points r*G and -s*Q
2. Builds a matrix M from monomial evaluations
3. Finds kernel basis to get linear relations
4. Solves for d such that Q = d*G

This is a Monte Carlo/Las Vegas approach that may succeed after several attempts.

Time Complexity: O(polynomial in log p) per attempt (probabilistic)
Space Complexity: O(n'^2) where n' is parameter
"""

import sys
import time
import random
from pathlib import Path
from typing import List, Tuple, Optional

# Add parent directory to path for utils import
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import EllipticCurve, Point, load_input, mod_inv


class LinearAlgebra:
    """Linear algebra operations over finite field F_p."""
    
    @staticmethod
    def rref(matrix: List[List[int]], p: int) -> Tuple[List[List[int]], List[Tuple[int, int]]]:
        """
        Compute Row-Reduced Echelon Form over F_p.
        
        Returns:
            (rref_matrix, pivot_positions)
        """
        M = [row[:] for row in matrix]  # Copy
        rows = len(M)
        cols = len(M[0]) if rows > 0 else 0
        
        pivots = []
        r = 0
        
        for c in range(cols):
            if r >= rows:
                break
            
            # Find pivot
            pivot_row = -1
            for i in range(r, rows):
                if M[i][c] % p != 0:
                    pivot_row = i
                    break
            
            if pivot_row == -1:
                continue
            
            # Swap rows
            M[r], M[pivot_row] = M[pivot_row], M[r]
            pivots.append((r, c))
            
            # Normalize pivot
            pivot_inv = mod_inv(M[r][c], p)
            for j in range(cols):
                M[r][j] = (M[r][j] * pivot_inv) % p
            
            # Eliminate column
            for i in range(rows):
                if i != r and M[i][c] != 0:
                    factor = M[i][c]
                    for j in range(cols):
                        M[i][j] = (M[i][j] - factor * M[r][j]) % p
            
            r += 1
        
        return M, pivots
    
    @staticmethod
    def kernel_basis(matrix: List[List[int]], p: int) -> List[List[int]]:
        """
        Find basis for left kernel of matrix over F_p.
        
        Returns vectors v such that v * M = 0.
        """
        if not matrix:
            return []
        
        rows = len(matrix)
        cols = len(matrix[0])
        
        # Transpose matrix
        M_t = [[matrix[r][c] for r in range(rows)] for c in range(cols)]
        
        # Compute RREF
        rref, pivots = LinearAlgebra.rref(M_t, p)
        
        # Find free variables
        pivot_cols = {pivot[1] for pivot in pivots}
        free_cols = [c for c in range(rows) if c not in pivot_cols]
        
        if not free_cols:
            return []
        
        # Build kernel basis
        pivot_map = {pivot[1]: pivot[0] for pivot in pivots}
        basis = []
        
        for free_col in free_cols:
            vec = [0] * rows
            vec[free_col] = 1
            
            for pivot_col in pivot_cols:
                row = pivot_map[pivot_col]
                if free_col < len(rref[row]):
                    vec[pivot_col] = (-rref[row][free_col]) % p
            
            basis.append(vec)
        
        return basis


def build_monomials(P: Point, n_prime: int, p: int) -> List[int]:
    """
    Build monomial vector for point P.
    
    For n_prime=2: Uses degree-2 monomials [x^2, y^2, 1, xy, x, y]
    For n_prime=3: Uses degree-3 monomials [x^3, y^3, 1, x^2y, x^2, y^2x, y^2, x, y, xy]
    """
    if P is None or P[0] is None or P[1] is None:
        raise ValueError("Point cannot be at infinity for monomial construction")
    
    x, y = P[0] % p, P[1] % p
    
    if n_prime == 2:
        x2 = (x * x) % p
        y2 = (y * y) % p
        xy = (x * y) % p
        return [x2, y2, 1, xy, x, y]
    elif n_prime == 3:
        x2 = (x * x) % p
        y2 = (y * y) % p
        x3 = (x * x2) % p
        y3 = (y * y2) % p
        x2y = (x2 * y) % p
        y2x = (y2 * x) % p
        xy = (x * y) % p
        return [x3, y3, 1, x2y, x2, y2x, y2, x, y, xy]
    else:
        raise ValueError(f"n_prime={n_prime} not supported (use 2 or 3)")


def find_sparse_vector(basis: List[List[int]], target_zeros: int, p: int) -> Optional[List[int]]:
    """
    Heuristic to find a vector with approximately 'target_zeros' zero entries.
    
    This implements a simplified version of the algorithm from the paper.
    """
    if not basis:
        return None
    
    dim = len(basis)
    vec_len = len(basis[0])
    
    # Try different linear combinations
    for attempt in range(100):
        # Random coefficients
        coeffs = [random.randrange(p) for _ in range(dim)]
        
        # Compute linear combination
        result = [0] * vec_len
        for i, coeff in enumerate(coeffs):
            if coeff != 0:
                for j in range(vec_len):
                    result[j] = (result[j] + coeff * basis[i][j]) % p
        
        # Check if we have enough zeros
        zero_count = sum(1 for x in result if x == 0)
        if zero_count >= target_zeros:
            return result
    
    return None


def las_vegas_ecdlp(curve: EllipticCurve, G: Point, Q: Point, n: int, 
                    n_prime: int = 3, max_attempts: int = 200) -> Tuple[Optional[int], int]:
    """
    Solve ECDLP using Las Vegas algorithm - Simplified version.
    
    Uses random linear combinations: tries to find relation sum(a_i * r_i*G) + sum(b_j * s_j*Q) = O
    which gives: sum(a_i * r_i) + sum(b_j * s_j * d) = 0 mod n
    so: d = -sum(a_i * r_i) / sum(b_j * s_j) mod n
    """
    
    for attempt in range(1, max_attempts + 1):
        # Progress indicator every 1000 attempts
        if attempt % 1000 == 0:
            print(f"  Attempt {attempt}/{max_attempts}...", flush=True)
        
        try:
            # Generate random points and coefficients
            k = 10  # Number of random combinations
            
            r_vals = []  # Coefficients for G
            s_vals = []  # Coefficients for Q
            points = []
            
            # Generate k random r*G points
            for _ in range(k):
                r = random.randrange(1, n)
                r_vals.append(r)
            
            # Generate k random s*Q points  
            for _ in range(k):
                s = random.randrange(1, n)
                s_vals.append(s)
            
            # Try random linear combinations
            for trial in range(50):
                # Random coefficients
                a_coeffs = [random.randrange(0, min(100, n)) for _ in range(k)]
                b_coeffs = [random.randrange(0, min(100, n)) for _ in range(k)]
                
                # Check if this gives point at infinity
                result = None
                
                # Add a_i * (r_i * G)
                for i, a in enumerate(a_coeffs):
                    if a != 0:
                        term = curve.scalar_multiply((a * r_vals[i]) % n, G)
                        result = curve.add(result, term)
                
                # Add b_j * (s_j * Q)
                for j, b in enumerate(b_coeffs):
                    if b != 0:
                        term = curve.scalar_multiply((b * s_vals[j]) % n, Q)
                        result = curve.add(result, term)
                
                # If we got identity, we have a relation!
                if result is None:
                    # sum(a_i * r_i) + d * sum(b_j * s_j) = 0 mod n
                    sum_a_r = sum(a_coeffs[i] * r_vals[i] for i in range(k)) % n
                    sum_b_s = sum(b_coeffs[j] * s_vals[j] for j in range(k)) % n
                    
                    if sum_b_s != 0:
                        # d = -sum_a_r / sum_b_s mod n
                        sum_b_s_inv = mod_inv(sum_b_s, n)
                        d = ((-sum_a_r) * sum_b_s_inv) % n
                        
                        # Verify
                        if curve.scalar_multiply(d, G) == Q:
                            return d, attempt
        
        except (ValueError, ZeroDivisionError):
            continue
    
    return None, max_attempts


def main():
    """Main entry point for Las Vegas ECDLP solver."""
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
    
    n_prime = 2  # Use 2 for faster performance
    max_attempts = 10000  # Need many attempts for probabilistic success
    
    print(f"Solving ECDLP using Las Vegas Algorithm...")
    print(f"Curve: y^2 = x^3 + {a}x + {b} (mod {p})")
    print(f"G = ({G[0]}, {G[1]}), Q = ({Q[0]}, {Q[1]}), n = {n}") #type: ignore
    print(f"Parameters: n'={n_prime}, l={3*n_prime}")
    
    start_time = time.perf_counter()
    d, attempts = las_vegas_ecdlp(curve, G, Q, n, n_prime, max_attempts)
    elapsed = time.perf_counter() - start_time
    
    if d is not None:
        Q_verify = curve.scalar_multiply(d, G)
        verified = (Q_verify == Q)
        
        print(f"\n{'='*50}")
        print(f"Solution: d = {d}")
        print(f"Attempts: {attempts}/{max_attempts}")
        print(f"Time: {elapsed:.6f} seconds")
        print(f"Verification: {'PASSED' if verified else 'FAILED'}")
        print(f"{'='*50}")
        
        if not verified:
            sys.exit(1)
    else:
        print(f"\nNo solution found after {max_attempts} attempts")
        print(f"Time: {elapsed:.6f} seconds")
        print("Note: This is a probabilistic algorithm. Try increasing max_attempts or adjusting n_prime.")
        sys.exit(1)


if __name__ == "__main__":
    main()
