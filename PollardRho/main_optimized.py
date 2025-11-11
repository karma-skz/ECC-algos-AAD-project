"""
Pollard Rho Algorithm for ECDLP - Optimized with C++

Solves the Elliptic Curve Discrete Logarithm Problem using a probabilistic
collision-finding approach with C++ backend for scalar multiplication.

Time Complexity: O(sqrt(n)) expected
Space Complexity: O(1)
"""

import sys
import time
import random
from pathlib import Path
from typing import List, Tuple, Optional
from math import gcd
import ctypes

# Add parent directory to path for utils import
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import EllipticCurve, Point, load_input, mod_inv


# Try to load C++ library
USE_CPP = False
ecc_lib = None

try:
    lib_path = Path(__file__).parent / 'ecc_fast.so'
    if lib_path.exists():
        ecc_lib = ctypes.CDLL(str(lib_path))
        
        # Configure function signatures
        ecc_lib.scalar_mult.argtypes = [
            ctypes.c_int64, ctypes.c_int64, ctypes.c_int64,  # k, gx, gy
            ctypes.c_int64, ctypes.c_int64, ctypes.c_int64,  # a, b, p
            ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)  # result_x, result_y
        ]
        ecc_lib.scalar_mult.restype = ctypes.c_int
        
        ecc_lib.point_add.argtypes = [
            ctypes.c_int64, ctypes.c_int64, ctypes.c_int64, ctypes.c_int64,  # x1, y1, x2, y2
            ctypes.c_int64, ctypes.c_int64,  # a, p
            ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)  # x3, y3
        ]
        ecc_lib.point_add.restype = ctypes.c_int
        
        USE_CPP = True
        print("✓ Using C++ optimized scalar multiplication")
except Exception as e:
    print(f"⚠ C++ library not available, using Python fallback: {e}")


def fast_scalar_mult(k: int, G: Point, curve: EllipticCurve) -> Point:
    """Fast scalar multiplication using C++ if available."""
    if not USE_CPP or G is None:
        return curve.scalar_multiply(k, G)
    
    gx, gy = G
    result_x = ctypes.c_int64()
    result_y = ctypes.c_int64()
    
    valid = ecc_lib.scalar_mult(
        ctypes.c_int64(k),
        ctypes.c_int64(gx),
        ctypes.c_int64(gy),
        ctypes.c_int64(curve.a),
        ctypes.c_int64(curve.b),
        ctypes.c_int64(curve.p),
        ctypes.byref(result_x),
        ctypes.byref(result_y)
    )
    
    if valid == 0:
        return None
    return (result_x.value, result_y.value)


def fast_point_add(P: Point, Q: Point, curve: EllipticCurve) -> Point:
    """Fast point addition using C++ if available."""
    if not USE_CPP:
        return curve.add(P, Q)
    
    if P is None:
        return Q
    if Q is None:
        return P
    
    x1, y1 = P
    x2, y2 = Q
    result_x = ctypes.c_int64()
    result_y = ctypes.c_int64()
    
    valid = ecc_lib.point_add(
        ctypes.c_int64(x1),
        ctypes.c_int64(y1),
        ctypes.c_int64(x2),
        ctypes.c_int64(y2),
        ctypes.c_int64(curve.a),
        ctypes.c_int64(curve.p),
        ctypes.byref(result_x),
        ctypes.byref(result_y)
    )
    
    if valid == 0:
        return None
    return (result_x.value, result_y.value)


def canonicalize(X: Point, A: int, B: int, curve: EllipticCurve, n: int) -> Tuple[Point, int, int]:
    """
    Apply negation map to choose canonical representative.
    Choose between X and -X based on y-coordinate to speed up collision detection.
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
    Each attempt gets a fresh random table for better coverage.
    """
    coefficients = []
    points = []
    
    for _ in range(m):
        # Random coefficients but avoid zero-zero
        u = random.randrange(1, n)  # At least 1
        v = random.randrange(0, n)
        
        # Compute R_i = u*G + v*Q
        uG = fast_scalar_mult(u, G, curve)
        vQ = fast_scalar_mult(v, Q, curve)
        R_i = fast_point_add(uG, vQ, curve)
        
        coefficients.append((u, v))
        points.append(R_i)
    
    return coefficients, points


def pollard_rho_optimized(curve: EllipticCurve, G: Point, Q: Point, n: int, 
                          max_steps: int = 10_000_000, partition_m: int = 32) -> Tuple[Optional[int], int]:
    """
    Solve ECDLP using Pollard's Rho with Floyd cycle detection and C++ optimization.
    Each call builds a fresh partition table for better randomness.
    """
    if G is None or Q is None:
        return None, 0
    
    # Build partition table (fresh random table each attempt)
    coefficients, points = make_partition_table(curve, G, Q, n, partition_m)
    
    def random_state() -> Tuple[Point, int, int]:
        """Generate random starting state."""
        A = random.randrange(n)
        B = random.randrange(n)
        X = fast_point_add(fast_scalar_mult(A, G, curve), fast_scalar_mult(B, Q, curve), curve)
        return canonicalize(X, A, B, curve, n)
    
    def step(state: Tuple[Point, int, int]) -> Tuple[Point, int, int]:
        """Perform one step of the random walk."""
        X, A, B = state
        
        # Choose partition based on x-coordinate
        idx = 0 if X is None else (X[0] % partition_m)
        R_i = points[idx]
        u_i, v_i = coefficients[idx]
        
        # Update: X' = X + R_i, A' = A + u_i, B' = B + v_i
        X_new = fast_point_add(X, R_i, curve)
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
        
        # Progress indicator every 50k steps
        if steps % 50000 == 0:
            print(f"  Step {steps}/{max_steps}...")
        
        X_t, A_t, B_t = tortoise
        X_h, A_h, B_h = hare
        
        if X_t == X_h:
            # Collision: A_t*G + B_t*Q = A_h*G + B_h*Q
            num = (A_t - A_h) % n
            den = (B_h - B_t) % n
            
            # Skip useless collisions
            if num == 0 and den == 0:
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
                    if fast_scalar_mult(d, G, curve) == Q:
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
                        
                        # Try all lifts
                        for k in range(g):
                            d_candidate = (d_base + k * n_reduced) % n
                            if fast_scalar_mult(d_candidate, G, curve) == Q:
                                return d_candidate, steps
                    except ZeroDivisionError:
                        pass
            
            # Failed collision - restart
            tortoise = random_state()
            hare = step(tortoise)
    
    return None, steps


def main():
    """Main entry point for Pollard Rho ECDLP solver."""
    script_dir = Path(__file__).parent
    
    if len(sys.argv) > 1:
        input_path = Path(sys.argv[1])
    else:
        input_path = script_dir / 'input' / 'small_test.txt'
    
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
    
    # Configuration - many short attempts work better than few long ones
    partition_m = 20  # Fewer partitions for simpler walks
    max_steps = 1_000_000  # Steps per attempt (1M)
    max_attempts = 20  # Many attempts to account for randomness
    
    print(f"Solving ECDLP using Optimized Pollard Rho...")
    print(f"Curve: y^2 = x^3 + {a}x + {b} (mod {p})")
    print(f"G = ({G[0]}, {G[1]}), Q = ({Q[0]}, {Q[1]}), n = {n}")
    print(f"Partitions: {partition_m}, Max steps/attempt: {max_steps}, Attempts: {max_attempts}")
    
    start_time = time.perf_counter()
    total_steps = 0
    d = None
    
    for attempt in range(1, max_attempts + 1):
        print(f"\n=== Attempt {attempt}/{max_attempts} ===")
        d_try, steps = pollard_rho_optimized(curve, G, Q, n, max_steps, partition_m)
        total_steps += steps
        
        if d_try is not None:
            d = d_try
            print(f"\n✓ Solution found on attempt {attempt}!")
            break
    
    elapsed = time.perf_counter() - start_time
    
    if d is not None:
        Q_verify = fast_scalar_mult(d, G, curve)
        verified = (Q_verify == Q)
        
        print(f"\n{'='*50}")
        print(f"Solution: d = {d} (found in attempt {attempt}/{max_attempts})")
        print(f"Total steps: {total_steps:,}")
        print(f"Time: {elapsed:.6f} seconds")
        print(f"Verification: {'PASSED' if verified else 'FAILED'}")
        print(f"{'='*50}")
        
        if not verified:
            sys.exit(1)
    else:
        print(f"\n{'='*50}")
        print(f"No solution found after {max_attempts} attempts")
        print(f"Total steps: {total_steps:,}")
        print(f"Time: {elapsed:.6f} seconds")
        print(f"Note: Pollard Rho is probabilistic and may need more attempts.")
        print(f"{'='*50}")
        sys.exit(1)


if __name__ == "__main__":
    main()
