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
                    n_prime: int = 2, max_attempts: int = 50) -> Tuple[Optional[int], int]:
    """
    Solve ECDLP using Las Vegas algorithm.
    
    Args:
        curve: The elliptic curve
        G: Base point (generator)
        Q: Target point (Q = d*G for unknown d)
        n: Order of base point G
        n_prime: Algorithm parameter (2 or 3)
        max_attempts: Maximum number of attempts
    
    Returns:
        (d, attempts) where d is the discrete log, or (None, attempts)
    """
    l = 3 * n_prime
    
    for attempt in range(1, max_attempts + 1):
        try:
            # Build matrix M from random point evaluations
            M = []
            I_list = []  # Multipliers for G
            J_list = []  # Multipliers for Q
            
            # Generate (3n'-1) rows from r*G
            for _ in range(3 * n_prime - 1):
                r = random.randrange(1, n)
                P = curve.scalar_multiply(r, G)
                if P is not None:
                    M.append(build_monomials(P, n_prime, curve.p))
                    I_list.append(r)
            
            # Generate (l+1) rows from -r*Q
            for _ in range(l + 1):
                r = random.randrange(1, n)
                P = curve.scalar_multiply(-r, Q)
                if P is not None:
                    M.append(build_monomials(P, n_prime, curve.p))
                    J_list.append(r)
            
            # Compute kernel basis
            kernel = LinearAlgebra.kernel_basis(M, curve.p)
            
            if not kernel:
                continue
            
            # Find sparse vector
            v = find_sparse_vector(kernel, l, curve.p)
            
            if v is None:
                continue
            
            # Extract solution
            A = 0
            B = 0
            
            for i in range(len(I_list)):
                if v[i] != 0:
                    A = (A + I_list[i]) % n
            
            for i in range(len(J_list)):
                v_idx = len(I_list) + i
                if v_idx < len(v) and v[v_idx] != 0:
                    B = (B + J_list[i]) % n
            
            # Compute d = A * B^(-1) mod n
            if B != 0:
                B_inv = mod_inv(B, n)
                d = (A * B_inv) % n
                
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
    
    n_prime = 2
    max_attempts = 100
    
    print(f"Solving ECDLP using Las Vegas Algorithm...")
    print(f"Curve: y^2 = x^3 + {a}x + {b} (mod {p})")
    print(f"G = ({G[0]}, {G[1]}), Q = ({Q[0]}, {Q[1]}), n = {n}")
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

# ==============================================================================
# 
#  Part 1: Modular Math Helpers
#
# ==============================================================================

def extended_gcd(a: int, b: int) -> Tuple[int, int, int]:
    """Returns (gcd, x, y) s.t. a*x + b*y = gcd"""
    if a == 0:
        return (b, 0, 1)
    else:
        gcd, x1, y1 = extended_gcd(b % a, a)
        x = y1 - (b // a) * x1
        y = x1
        return (gcd, x, y)

def mod_inv(a: int, m: int) -> int:
    """Calculates the modular inverse of a (mod m)"""
    gcd, x, y = extended_gcd(a, m)
    if gcd != 1:
        raise ValueError(f'Modular inverse does not exist for {a} mod {m}')
    else:
        return x % m

# ==============================================================================
#
#  Part 2: Elliptic Curve and Point Classes (Affine Coords)
#
# ==============================================================================

class Point:
    """A point on an elliptic curve [x, y] or the identity O (None)"""
    
    def __init__(self, x: Optional[int], y: Optional[int]):
        self.x = x
        self.y = y

    def is_identity(self) -> bool:
        """Check if the point is the identity"""
        return self.x is None or self.y is None

    def __eq__(self, other) -> bool:
        if not isinstance(other, Point):
            return False
        return self.x == other.x and self.y == other.y

    def __str__(self) -> str:
        if self.is_identity():
            return "Point(O)"
        return f"Point({self.x}, {self.y})"
        
    def __hash__(self):
        return hash((self.x, self.y))

class EllipticCurve:
    """
    Represents an elliptic curve y^2 = x^3 + ax + b (mod p)
    
    Args:
        a, b: Curve coefficients
        p: The prime modulus of the finite field
    """
    
    def __init__(self, a: int, b: int, p: int):
        self.a = a
        self.b = b
        self.p = p
        # Validate curve
        if (4 * a**3 + 27 * b**2) % p == 0:
            raise ValueError("Curve is singular")
    def is_on_curve(self, P: 'Point') -> bool:
        if P.is_identity():
            return True
        x, y = P.x % self.p, P.y % self.p # type: ignore
        return (y * y - (x * x * x + self.a * x + self.b)) % self.p == 0

    def identity(self) -> Point:
        """Returns the identity point O"""
        return Point(None, None)

    def negate(self, p1: Point) -> Point:
        """Returns -P"""
        if p1.is_identity():
            return self.identity()
        return Point(p1.x, -p1.y % self.p) # type: ignore

    def add(self, p1: Point, p2: Point) -> Point:
        """Adds two points P1 + P2 on the curve"""
        if p1.is_identity():
            return p2
        if p2.is_identity():
            return p1

        # Handle P1 == P2 (doubling)
        if p1 == p2:
            return self.double(p1)

        # Handle P1 == -P2
        if p1.x == p2.x and p1.y != p2.y:
            return self.identity()

        # Standard addition
        x1, y1 = p1.x, p1.y
        x2, y2 = p2.x, p2.y

        # Calculate slope 's'
        try:
            # s = (y2 - y1) * mod_inv(x2 - x1, p)
            s_num = (y2 - y1) % self.p # type: ignore
            s_den = mod_inv((x2 - x1) % self.p, self.p) # type: ignore
            s = (s_num * s_den) % self.p
        except ValueError:
            # This case (x1=x2, y1!=y2) is already handled,
            # but a duplicate x1=x2,y1=y2 should go to double()
            return self.identity()

        # Calculate new point
        x3 = (s**2 - x1 - x2) % self.p # type: ignore
        y3 = (s * (x1 - x3) - y1) % self.p # type: ignore

        return Point(x3, y3)

    def double(self, p1: Point) -> Point:
        """Doubles a point P (P + P)"""
        if p1.is_identity() or p1.y == 0:
            return self.identity()
            
        x, y = p1.x, p1.y
        
        # Calculate slope 's'
        # s = (3*x^2 + a) * mod_inv(2*y, p)
        s_num = (3 * x**2 + self.a) % self.p # type: ignore
        s_den = mod_inv((2 * y) % self.p, self.p) # type: ignore
        s = (s_num * s_den) % self.p

        # Calculate new point
        x3 = (s**2 - 2 * x) % self.p # type: ignore
        y3 = (s * (x - x3) - y) % self.p # type: ignore
        
        return Point(x3, y3)

    def multiply(self, p_in: Point, n: int) -> Point:
        """Scalar multiplication n*P using double-and-add"""
        result = self.identity()
        current = p_in
        
        if n < 0:
            n = abs(n)
            current = self.negate(current)

        while n > 0:
            if n % 2 == 1:
                result = self.add(result, current)
            current = self.double(current)
            n //= 2
            
        return result

# ==============================================================================
#
#  Part 3: Linear Algebra over a Finite Field (F_p)
#
# ==============================================================================

class LinearAlgebraFp:
    """Contains static methods for linear algebra over F_p"""

    @staticmethod
    def _create_matrix_copy(M: List[List[int]]) -> List[List[int]]:
        return [[cell for cell in row] for row in M]

    @staticmethod
    def find_left_kernel_basis(M: List[List[int]], p: int) -> List[List[int]]:
        """
        Computes the basis for the left kernel of M over F_p.
        This is the standard (right) kernel of M_transpose.
        """
        if not M:
            return []
        
        num_rows_M = len(M)
        num_cols_M = len(M[0])
        
        # Transpose M
        Mt = [[M[r][c] for r in range(num_rows_M)] for c in range(num_cols_M)]
        
        # Compute RREF of M_transpose
        rref, pivots = LinearAlgebraFp.compute_rref(Mt, p)
        
        # Extract kernel basis
        kernel_basis = []
        pivot_cols = [p[1] for p in pivots]
        free_var_cols = [c for c in range(num_rows_M) if c not in pivot_cols]
        
        pivot_row_map = {p[1]: p[0] for p in pivots}

        for free_col in free_var_cols:
            vec = [0] * num_rows_M
            vec[free_col] = 1
            for pivot_col in pivot_cols:
                row = pivot_row_map[pivot_col]
                if row < len(rref) and free_col < len(rref[row]):
                    vec[pivot_col] = -rref[row][free_col] % p
            kernel_basis.append(vec)
            
        return kernel_basis

    @staticmethod
    def compute_rref(A: List[List[int]], p: int) -> Tuple[List[List[int]], List[Tuple[int, int]]]:
        """Computes the Row-Reduced Echelon Form (RREF) of matrix A over F_p"""
        M = LinearAlgebraFp._create_matrix_copy(A)
        rows = len(M)
        cols = len(M[0]) if rows > 0 else 0
        
        pivots = []
        r = 0
        c = 0
        
        while r < rows and c < cols:
            # Find pivot
            i_pivot = -1
            for i in range(r, rows):
                if M[i][c] != 0:
                    i_pivot = i
                    break
            
            if i_pivot == -1:
                c += 1
                continue

            # Swap rows
            M[r], M[i_pivot] = M[i_pivot], M[r]
            pivots.append((r, c))

            # Normalize pivot row
            pivot_val = M[r][c]
            inv_pivot = mod_inv(pivot_val, p)
            for j in range(c, cols):
                M[r][j] = (M[r][j] * inv_pivot) % p
                
            # Eliminate other rows
            for i in range(rows):
                if i == r:
                    continue
                factor = M[i][c]
                if factor == 0:
                    continue
                for j in range(c, cols):
                    M[i][j] = (M[i][j] - factor * M[r][j]) % p
            
            r += 1
            c += 1
            
        return M, pivots

    @staticmethod
    def solve_problem_l_heuristic(K_basis: List[List[int]], l: int, p: int) -> Optional[List[int]]:
        """
        Implements Algorithm 2 from the paper.
        This is a *heuristic* to find a vector with 'l' zeros.
        """
        if not K_basis:
            return None
        
        # Total columns = 2l (since 3n' = l)
        total_cols = len(K_basis[0])
        
        # --- Attempt 1: Diagonalize first block (cols 0 to l-1) ---
        K_copy1 = LinearAlgebraFp._create_matrix_copy(K_basis)
        
        # Perform Gaussian elimination on the first block
        # This is a partial RREF, not a full one
        r = 0
        for c in range(l):
            if r >= l: break
            
            # Find pivot
            i_pivot = -1
            for i in range(r, l):
                if K_copy1[i][c] != 0:
                    i_pivot = i
                    break
            if i_pivot == -1: continue

            # Swap
            K_copy1[r], K_copy1[i_pivot] = K_copy1[i_pivot], K_copy1[r]
            
            # Normalize
            inv_pivot = mod_inv(K_copy1[r][c], p)
            for j in range(total_cols):
                K_copy1[r][j] = (K_copy1[r][j] * inv_pivot) % p
                
            # Eliminate
            for i in range(l):
                if i == r: continue
                factor = K_copy1[i][c]
                for j in range(total_cols):
                    K_copy1[i][j] = (K_copy1[i][j] - factor * K_copy1[r][j]) % p
            r += 1
        
        # Check for a solution
        for row in K_copy1:
            if row[l:total_cols].count(0) > 0:
                # This row has at least l-1 zeros from the diagonal block
                # and at least 1 zero in the second block.
                # Let's count total zeros to be sure.
                if row.count(0) >= l:
                    return row

        # --- Attempt 2: Diagonalize second block (cols l to 2l-1) ---
        K_copy2 = LinearAlgebraFp._create_matrix_copy(K_basis)
        
        r = 0
        for c in range(l, total_cols):
            if r >= l: break
            
            i_pivot = -1
            for i in range(r, l):
                if K_copy2[i][c] != 0:
                    i_pivot = i
                    break
            if i_pivot == -1: continue

            K_copy2[r], K_copy2[i_pivot] = K_copy2[i_pivot], K_copy2[r]
            
            inv_pivot = mod_inv(K_copy2[r][c], p)
            for j in range(total_cols):
                K_copy2[r][j] = (K_copy2[r][j] * inv_pivot) % p
                
            for i in range(l):
                if i == r: continue
                factor = K_copy2[i][c]
                for j in range(total_cols):
                    K_copy2[i][j] = (K_copy2[i][j] - factor * K_copy2[r][j]) % p
            r += 1

        # Check for a solution
        for row in K_copy2:
            if row[0:l].count(0) > 0:
                if row.count(0) >= l:
                    return row
                    
        # No solution found by this heuristic
        return None


# ==============================================================================
#
#  Part 4: The Las Vegas ECDLP Algorithm
#
# ==============================================================================

class LasVegasECDLP:
    """
    Implements the full algorithm based on the paper.
    
    Args:
        curve: An EllipticCurve object
        G: The base point (generator)
        Q: The public point (Q = d*G)
        n: The order of the subgroup generated by G
    """
    
    def __init__(self, curve: EllipticCurve, G: Point, Q: Point, n: int):
        self.curve = curve
        self.G = G
        self.Q = Q
        self.n = n  # Group order
        self.p = curve.p # Field prime
        self.rng = random.SystemRandom()

    def _build_monomial_row(self, point: Point, n_prime: int) -> List[int]:
        """
        [cite_start]Creates a row for matrix M from a point [cite: 75, 79]
        Assumes projective coords (x, y, 1)
        """
        x, y, z = point.x, point.y, 1
        
        if n_prime == 2:
            # Monomials: x^2, y^2, z^2, xy, xz, yz
            # With z=1: x^2, y^2, 1, xy, x, y
            return [
                (x*x) % self.p, # type: ignore
                (y*y) % self.p, # type: ignore
                1,
                (x*y) % self.p, # type: ignore
                x,
                y
            ]
        elif n_prime == 3:
            # Monomials: x^3, y^3, z^3, x^2y, x^2z, y^2x, y^2z, z^2x, z^2y, xyz
            # With z=1: x^3, y^3, 1, x^2y, x^2, y^2x, y^2, x, y, xy
            x2, y2 = (x*x) % self.p, (y*y) % self.p # type: ignore
            x3, y3 = (x*x2) % self.p, (y*y2) % self.p # type: ignore
            return [
                x3, y3, 1, # type: ignore
                (x2*y) % self.p, x2, # type: ignore
                (y2*x) % self.p, y2, # type: ignore
                x, y, 
                (x*y) % self.p # type: ignore
            ]
        else:
            # This can be generalized, but hard-coded for now
            raise ValueError(f"n_prime = {n_prime} not supported. Use 2 or 3.")

    def _build_matrix_M(self, n_prime: int, l: int) -> Tuple[List[List[int]], List[int], List[int]]:
        """
        [cite_start]Generates random points and builds the matrix M [cite: 133-151]
        """
        num_cols = (n_prime + 1) * (n_prime + 2) // 2
        num_rows_G = 3 * n_prime - 1
        num_rows_Q = l + 1
        
        M = []
        I_list = [] # List of r for r*G
        J_list = [] # List of r for -r*Q
        
        used_r_G = set()
        used_r_Q = set()

        # Generate (3n' - 1) rows from r*G
        for _ in range(num_rows_G):
            point_rG = self.curve.identity()
            r = -1 # Init to invalid value
            
            # Loop until we get a non-identity point
            while point_rG.is_identity():
                r = self.rng.randrange(1, self.n)
                while r in used_r_G:
                    r = self.rng.randrange(1, self.n)
                
                point_rG = self.curve.multiply(self.G, r)
            
            # We now have a valid r and a valid (non-identity) point
            used_r_G.add(r)
            I_list.append(r)
            M.append(self._build_monomial_row(point_rG, n_prime))


        # Generate (l + 1) rows from -r*Q
        for _ in range(num_rows_Q):
            point_minus_rQ = self.curve.identity()
            r = -1 # Init to invalid value

            # Loop until we get a non-identity point
            while point_minus_rQ.is_identity():
                r = self.rng.randrange(1, self.n)
                while r in used_r_Q:
                    r = self.rng.randrange(1, self.n)
                
                point_minus_rQ = self.curve.multiply(self.Q, -r)

            # We now have a valid r and a valid (non-identity) point
            used_r_Q.add(r)
            J_list.append(r)
            M.append(self._build_monomial_row(point_minus_rQ, n_prime))
            
        return M, I_list, J_list

    def solve(self, n_prime: int, max_tries: int = 50) -> Tuple[Optional[int], int]:
        """
        Attempts to solve the ECDLP
        
        Args:
            n_prime: The 'n'' parameter from the paper. Small integer (e.g., 2 or 3)
            max_tries: Number of times to repeat Algorithm 1
            
        Returns:
            A tuple (found_d, attempts):
            - found_d: The discovered secret (int) or None if not found.
            - attempts: The number of attempts taken.
        """
        
        l = 3 * n_prime
        # Keep noisy prints minimal; high-level metadata printed by main
        
        for attempt in range(1, max_tries + 1):
            # 1. [cite_start]Build Matrix M and get random multipliers [cite: 133-151]
            M, I_list, J_list = self._build_matrix_M(n_prime, l)
            
            # 2. [cite_start]Compute Left-Kernel K of M [cite: 152]
            K_basis = LinearAlgebraFp.find_left_kernel_basis(M, self.p)
            
            if not K_basis:
                continue
                
            # 3. [cite_start]Solve Problem L (Heuristic) [cite: 153]
            v = LinearAlgebraFp.solve_problem_l_heuristic(K_basis, l, self.p)
            
            if v is None:
                continue

            # 4. Solution Found! [cite_start]Calculate A and B [cite: 154-163]
            A = 0
            B = 0
            
            # Sum for G points (indices 0 to 3n'-2)
            for i in range(len(I_list)):
                if v[i] != 0:
                    A = (A + I_list[i]) % self.n
                    
            # Sum for Q points (indices 3n'-1 to 3n'+l)
            for i in range(len(J_list)):
                v_index = len(I_list) + i
                if v[v_index] != 0:
                    B = (B + J_list[i]) % self.n
            
            # 5. [cite_start]Final Calculation: d = A * B^-1 (mod n) [cite: 164]
            # Note: We use mod 'n' (group order) as per ECDLP
            # [cite_start]The paper uses 'p' [cite: 164] [cite_start]because it *defined* p as the group order [cite: 32]
            try:
                B_inv = mod_inv(B, self.n)
                d = (A * B_inv) % self.n
                
                # Verification
                Q_check = self.curve.multiply(self.G, d)
                if Q_check == self.Q:
                    return d, attempt # <-- CHANGED: Return d and attempt count
            
            except ValueError:
                pass
                
        return None, max_tries # <-- CHANGED: Return None and max attempts

# ==============================================================================
#
#  Part 5: Main Execution (Example)
#
# ==============================================================================


def load_positional_input(file_path: Path):
    with file_path.open('r') as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    if len(lines) < 5:
        raise ValueError("Input file must contain 5 non-empty lines: p | a b | Gx Gy | n | Qx Qy")
    def ints(s: str):
        return list(map(int, s.split()))
    p = ints(lines[0])[0]
    a, b = ints(lines[1])
    Gx, Gy = ints(lines[2])
    n = ints(lines[3])[0]
    Qx, Qy = ints(lines[4])
    return p, a, b, (Gx, Gy), n, (Qx, Qy)

if __name__ == "__main__":
    script_dir = Path(__file__).parent
    default_path = script_dir / 'input' / 'filename.txt'
    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else default_path
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    p, a, b, G_tuple, n, Q_tuple = load_positional_input(input_path)

    # Build curve and points
    curve = EllipticCurve(a, b, p)
    G = Point(G_tuple[0], G_tuple[1])
    Q = Point(Q_tuple[0], Q_tuple[1])

    # Sanity checks
    if not curve.is_on_curve(G):
        raise ValueError("G is not on the curve")
    if not curve.is_on_curve(Q):
        raise ValueError("Q is not on the curve")
    nG = curve.multiply(G, n)
    if not nG.is_identity():
        print("Warning: n*G != O; provided n may not be exact order")

    method = "LasVegas"
    N_PRIME = 2
    MAX_TRIES = 100

    start = time.perf_counter()
    solver = LasVegasECDLP(curve=curve, G=G, Q=Q, n=n)
    found_d, attempts = solver.solve(n_prime=N_PRIME, max_tries=MAX_TRIES)
    elapsed = time.perf_counter() - start

    # Output metadata
    print(f"Curve: p={p}, a={a}, b={b}")
    print(f"G=({G.x},{G.y}), Q=({Q.x},{Q.y}), n={n}")
    print(f"Method: {method}")
    print(f"Params: n'={N_PRIME}, l={3*N_PRIME}")
    print(f"Time elapsed: {elapsed:.6f} s")

    if found_d is not None:
        print(f"Found d: {found_d}")
        check = curve.multiply(G, found_d % n)
        match = (check == Q)
        print(f"Verified: {match}")
        print(f"Match: {match}")
    else:
        print("No solution found (None)")