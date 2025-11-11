"""
Optimized Las Vegas with C++ backend for scalar multiplication.

This uses ctypes to call compiled C++ code for faster elliptic curve operations.
To compile: g++ -O3 -shared -fPIC ecc_fast.cpp -o ecc_fast.so
"""

import sys
import time
import random
import ctypes
from pathlib import Path
from typing import Optional, Tuple

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import EllipticCurve, Point, load_input, mod_inv


# Try to load C++ library
try:
    lib_path = Path(__file__).parent / "ecc_fast.so"
    if lib_path.exists():
        ecc_lib = ctypes.CDLL(str(lib_path))
        
        # Define function signatures
        ecc_lib.scalar_mult.argtypes = [
            ctypes.c_longlong,  # k
            ctypes.c_longlong,  # Gx
            ctypes.c_longlong,  # Gy
            ctypes.c_longlong,  # a
            ctypes.c_longlong,  # b
            ctypes.c_longlong,  # p
            ctypes.POINTER(ctypes.c_longlong),  # result_x
            ctypes.POINTER(ctypes.c_longlong),  # result_y
        ]
        ecc_lib.scalar_mult.restype = ctypes.c_int  # 0 if infinity, 1 if valid point
        
        USE_CPP = True
        print("✓ Using C++ optimized scalar multiplication")
    else:
        USE_CPP = False
        print("⚠ C++ library not found, using Python fallback (slower)")
except Exception as e:
    USE_CPP = False
    print(f"⚠ C++ library load failed: {e}, using Python fallback")


def fast_scalar_mult(k: int, G: Point, curve: EllipticCurve) -> Optional[Point]:
    """Fast scalar multiplication using C++ if available."""
    if not USE_CPP or G is None:
        return curve.scalar_multiply(k, G)
    
    result_x = ctypes.c_longlong()
    result_y = ctypes.c_longlong()
    
    is_valid = ecc_lib.scalar_mult(
        k % curve.p,
        G[0],
        G[1],
        curve.a,
        curve.b,
        curve.p,
        ctypes.byref(result_x),
        ctypes.byref(result_y)
    )
    
    if is_valid:
        return (result_x.value, result_y.value)
    else:
        return None


def las_vegas_optimized(curve: EllipticCurve, G: Point, Q: Point, n: int, 
                       max_attempts: int = 100000) -> Tuple[Optional[int], int]:
    """
    Optimized Las Vegas using C++ for scalar multiplications.
    """
    
    for attempt in range(1, max_attempts + 1):
        if attempt % 1000 == 0:  # Progress every 1000 attempts
            print(f"  Attempt {attempt}/{max_attempts}...", flush=True)
        
        try:
            k = 10  # Number of random points
            
            r_vals = [random.randrange(1, n) for _ in range(k)]
            s_vals = [random.randrange(1, n) for _ in range(k)]
            
            # Try random linear combinations
            for trial in range(50):
                a_coeffs = [random.randrange(0, min(1000, n)) for _ in range(k)]
                b_coeffs = [random.randrange(0, min(1000, n)) for _ in range(k)]
                
                # Compute sum using fast multiplication
                result = None
                
                # Add a_i * (r_i * G)
                for i, a in enumerate(a_coeffs):
                    if a != 0:
                        term = fast_scalar_mult((a * r_vals[i]) % n, G, curve)
                        result = curve.add(result, term)
                
                # Add b_j * (s_j * Q)
                for j, b in enumerate(b_coeffs):
                    if b != 0:
                        term = fast_scalar_mult((b * s_vals[j]) % n, Q, curve)
                        result = curve.add(result, term)
                
                # Check if we hit identity
                if result is None:
                    sum_a_r = sum(a_coeffs[i] * r_vals[i] for i in range(k)) % n
                    sum_b_s = sum(b_coeffs[j] * s_vals[j] for j in range(k)) % n
                    
                    if sum_b_s != 0:
                        sum_b_s_inv = mod_inv(sum_b_s, n)
                        d = ((-sum_a_r) * sum_b_s_inv) % n
                        
                        # Verify with fast multiplication
                        if fast_scalar_mult(d, G, curve) == Q:
                            return d, attempt
        
        except (ValueError, ZeroDivisionError):
            continue
    
    return None, max_attempts


def main():
    """Main entry point."""
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
    
    print(f"Solving ECDLP using Optimized Las Vegas...")
    print(f"Curve: y^2 = x^3 + {a}x + {b} (mod {p})")
    print(f"G = ({G[0]}, {G[1]}), Q = ({Q[0]}, {Q[1]}), n = {n}")
    print()
    
    max_attempts = 50000  # Increased attempts
    
    start_time = time.perf_counter()
    d, attempts = las_vegas_optimized(curve, G, Q, n, max_attempts)
    elapsed = time.perf_counter() - start_time
    
    if d is not None:
        print(f"\n{'='*50}")
        print(f"Solution: d = {d}")
        print(f"Attempts: {attempts}/{max_attempts}")
        print(f"Time: {elapsed:.6f} seconds")
        
        # Verify
        Q_verify = fast_scalar_mult(d, G, curve)
        verified = (Q_verify == Q)
        print(f"Verification: {'PASSED' if verified else 'FAILED'}")
        print(f"{'='*50}")
        
        if not verified:
            sys.exit(1)
    else:
        print(f"\nNo solution found after {max_attempts} attempts")
        print(f"Time: {elapsed:.6f} seconds")
        print("Note: Las Vegas is highly probabilistic and may not succeed on large inputs.")
        sys.exit(1)


if __name__ == "__main__":
    main()
