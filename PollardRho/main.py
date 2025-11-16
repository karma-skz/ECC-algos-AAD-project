"""
Pollard Rho Algorithm for ECDLP

Solves the Elliptic Curve Discrete Logarithm Problem using a probabilistic
collision-finding approach. Given Q = d*G, finds d by:
1. Using pseudo-random walks on the curve with partition function
2. Applying Floyd's cycle detection to find collision
3. Using negation map to speed up collision detection

Time Complexity: O(sqrt(n)) expected
Space Complexity: O(1)
"""

import sys
import time
import random
from pathlib import Path
from typing import List, Tuple, Optional
from math import gcd

# Add parent directory to path for utils import
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import EllipticCurve, Point, load_input, mod_inv


def canonicalize(X: Point, A: int, B: int, curve: EllipticCurve, n: int) -> Tuple[Point, int, int]:
    """
    Apply negation map to choose canonical representative.
    
    Choose between X and -X based on y-coordinate to speed up collision detection.
    If we negate the point, we also negate the coefficients: -(A*G + B*Q) = (-A)*G + (-B)*Q
    
    Args:
        X: Point on curve
        A, B: Coefficients such that X = A*G + B*Q
        curve: The elliptic curve
        n: Order of base point
    
    Returns:
        Canonicalized (X, A, B)
    """
    if X is None:
        return X, A % n, B % n
    
    x, y = X
    # Choose representative with smaller y-coordinate
    if y > (curve.p - y) % curve.p:
        X = (x % curve.p, (-y) % curve.p)
        A = (-A) % n
        B = (-B) % n
    
    return X, A % n, B % n


def make_partition_table(curve: EllipticCurve, G: Point, Q: Point, n: int, m: int = 16) -> Tuple[List[Tuple[int, int]], List[Point]]:
    """
    Build m-way partition jump table for pseudo-random walks.
    
    For each partition i, choose random (u_i, v_i) and compute R_i = u_i*G + v_i*Q.
    On encountering partition i, we update: X <- X + R_i, A <- A + u_i, B <- B + v_i
    
    Args:
        curve: The elliptic curve
        G: Base point
        Q: Target point
        n: Order of base point
        m: Number of partitions
    
    Returns:
        (coefficients, points) where coefficients[i] = (u_i, v_i) and points[i] = R_i
    """
    coefficients = []
    points = []
    
    for _ in range(m):
        # Choose random non-zero coefficients
        while True:
            u = random.randrange(n)
            v = random.randrange(n)
            if (u | v) != 0:  # At least one non-zero
                break
        
        # Compute R_i = u*G + v*Q
        uG = curve.scalar_multiply(u, G)
        vQ = curve.scalar_multiply(v, Q)
        R_i = curve.add(uG, vQ)
        
        coefficients.append((u, v))
        points.append(R_i)
    
    return coefficients, points


def pollard_rho_ecdlp(curve: EllipticCurve, G: Point, Q: Point, n: int, 
                      max_steps: int = 2_000_000, partition_m: int = 32) -> Tuple[Optional[int], int]:
    """
    Solve ECDLP using Pollard's Rho with Floyd cycle detection.
    
    Args:
        curve: The elliptic curve
        G: Base point (generator)
        Q: Target point (Q = d*G for unknown d)
        n: Order of base point G
        max_steps: Maximum iterations before giving up
        partition_m: Number of partitions for random walk
    
    Returns:
        (d, steps) where d is the discrete log, or (None, steps) if not found
    """
    if G is None or Q is None:
        return None, 0
    
    # Build partition table
    coefficients, points = make_partition_table(curve, G, Q, n, partition_m)
    
    def random_state() -> Tuple[Point, int, int]:
        """Generate random starting state."""
        A = random.randrange(n)
        B = random.randrange(n)
        X = curve.add(curve.scalar_multiply(A, G), curve.scalar_multiply(B, Q))
        return canonicalize(X, A, B, curve, n)
    
    def step(state: Tuple[Point, int, int]) -> Tuple[Point, int, int]:
        """Perform one step of the random walk."""
        X, A, B = state
        
        # Choose partition based on x-coordinate
        idx = 0 if X is None else (X[0] % partition_m)
        R_i = points[idx]
        u_i, v_i = coefficients[idx]
        
        # Update: X' = X + R_i, A' = A + u_i, B' = B + v_i
        X_new = curve.add(X, R_i)
        A_new = (A + u_i) % n
        B_new = (B + v_i) % n
        
        return canonicalize(X_new, A_new, B_new, curve, n)
    
    # Floyd's cycle detection
    tortoise = random_state()
    hare = step(tortoise)
    steps = 0
    
    while steps < max_steps:
        tortoise = step(tortoise)
        hare = step(step(hare))
        steps += 1
        
        X_t, A_t, B_t = tortoise
        X_h, A_h, B_h = hare
        
        if X_t == X_h:
            # Collision: A_t*G + B_t*Q = A_h*G + B_h*Q
            # => (A_t - A_h)*G = (B_h - B_t)*Q = (B_h - B_t)*d*G
            # => A_t - A_h = d*(B_h - B_t) (mod n)
            
            num = (A_t - A_h) % n
            den = (B_h - B_t) % n
            
            # Skip useless collisions where both coefficients are equal
            if num == 0 and den == 0:
                # Restart with new random state
                tortoise = random_state()
                hare = step(tortoise)
                continue
            
            g = gcd(den, n)
            
            if g == 1 and den != 0:
                # Simple case: den is coprime to n
                try:
                    den_inv = mod_inv(den, n)
                    d = (num * den_inv) % n
                    
                    # Verify solution
                    if curve.scalar_multiply(d, G) == Q:
                        return d, steps
                except ZeroDivisionError:
                    pass
            elif g > 1 and num % g == 0:
                # Need to solve with gcd(den, n) > 1
                n_reduced = n // g
                den_reduced = (den // g) % n_reduced
                num_reduced = (num // g) % n_reduced
                
                if den_reduced != 0:
                    try:
                        den_inv = mod_inv(den_reduced, n_reduced)
                        d_base = (num_reduced * den_inv) % n_reduced
                        
                        # Try all lifts: d = d_base + k*(n/g) for k = 0, ..., g-1
                        for k in range(g):
                            d_candidate = (d_base + k * n_reduced) % n
                            if curve.scalar_multiply(d_candidate, G) == Q:
                                return d_candidate, steps
                    except ZeroDivisionError:
                        pass
            
            # Failed collision - restart with new random state
            tortoise = random_state()
            hare = step(tortoise)
    
    return None, steps


def main():
    """Main entry point for Pollard Rho ECDLP solver."""
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
    
    # Configuration
    partition_m = 32
    max_steps = 10_000_000  # Increased for better success rate
    max_attempts = 10
    
    print(f"Solving ECDLP using Pollard Rho...")
    print(f"Curve: y^2 = x^3 + {a}x + {b} (mod {p})")
    print(f"G = ({G[0]}, {G[1]}), Q = ({Q[0]}, {Q[1]}), n = {n}") # type: ignore
    print(f"Partitions: {partition_m}, Max steps/attempt: {max_steps}")
    
    start_time = time.perf_counter()
    total_steps = 0
    d = None
    
    for attempt in range(1, max_attempts + 1):
        d_try, steps = pollard_rho_ecdlp(curve, G, Q, n, max_steps, partition_m)
        total_steps += steps
        
        if d_try is not None:
            d = d_try
            print(f"Solution found on attempt {attempt}")
            break
    
    elapsed = time.perf_counter() - start_time
    
    if d is not None:
        Q_verify = curve.scalar_multiply(d, G)
        verified = (Q_verify == Q)
        
        print(f"\n{'='*50}")
        print(f"Solution: d = {d}")
        print(f"Total steps: {total_steps}")
        print(f"Time: {elapsed:.6f} seconds")
        print(f"Verification: {'PASSED' if verified else 'FAILED'}")
        print(f"{'='*50}")
        
        if not verified:
            sys.exit(1)
    else:
        print(f"\nNo solution found after {max_attempts} attempts")
        print(f"Total steps: {total_steps}")
        print(f"Time: {elapsed:.6f} seconds")
        sys.exit(1)


if __name__ == "__main__":
    main()
