"""
MOV Attack - Optimized with C++

Optimized version using C++ for scalar multiplication operations.
The pairing computation itself remains in Python as it's algorithm-specific.
"""

import sys
import time
import ctypes
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import EllipticCurve, Point, load_input
from MOV.main import (
    find_embedding_degree, generate_irreducible_poly,
    Polynomial, FieldElement, bsgs_field
)

# Load C++ library for fast scalar multiplication
USE_CPP = False
ecc_lib = None

try:
    lib_path = Path(__file__).parent.parent / 'utils' / 'cpp' / 'ecc_fast.so'
    if not lib_path.exists():
        lib_path = Path(__file__).parent / 'ecc_fast.so'
    
    if lib_path.exists():
        ecc_lib = ctypes.CDLL(str(lib_path))
        
        # Configure function signatures
        ecc_lib.scalar_mult.argtypes = [
            ctypes.c_int64, ctypes.c_int64, ctypes.c_int64,  # k, gx, gy
            ctypes.c_int64, ctypes.c_int64, ctypes.c_int64,  # a, b, p
            ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)  # result
        ]
        ecc_lib.scalar_mult.restype = ctypes.c_int
        
        USE_CPP = True
        print("✓ Using C++ optimized scalar multiplication")
except Exception as e:
    print(f"⚠ C++ library not available, using Python: {e}")


def fast_scalar_mult(k: int, G: Point, curve: EllipticCurve) -> Point:
    """Fast scalar multiplication using C++ if available."""
    if not USE_CPP or G is None:
        return curve.scalar_multiply(k, G)
    
    gx, gy = G
    result_x = ctypes.c_int64()
    result_y = ctypes.c_int64()
    
    valid = ecc_lib.scalar_mult( # type: ignore
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


def mov_attack_optimized(curve: EllipticCurve, G: Point, Q: Point, n: int, 
                         max_k: int = 12):
    """
    Optimized MOV attack using C++ for base operations.
    Pairing computation remains placeholder.
    """
    p = curve.p
    
    # Step 1: Find embedding degree
    print(f"Step 1: Finding embedding degree (max k={max_k})...")
    k = find_embedding_degree(p, n, max_k)
    
    if k is None:
        print(f"✗ Embedding degree > {max_k}, MOV attack not practical")
        return None
    
    print(f"✓ Found embedding degree: k = {k}")
    print(f"  Extension field size: p^{k} = {p**k}")
    
    # Step 2: Generate irreducible polynomial
    print(f"\nStep 2: Generating irreducible polynomial of degree {k}...")
    mod_poly = generate_irreducible_poly(k, p)
    print(f"✓ Irreducible polynomial generated")
    
    # Verify base points using optimized operations
    print(f"\nVerifying points using {'C++' if USE_CPP else 'Python'} operations...")
    nG = fast_scalar_mult(n, G, curve)
    if nG is not None:
        print("⚠ Warning: n*G ≠ O")
    
    # Step 3-5: Pairing computation (placeholder)
    print(f"\nStep 3-5: Computing pairings and solving DLP...")
    print("⚠ Note: Full pairing computation requires:")
    print("  • Extension field arithmetic for E(F_p^k)")
    print("  • Miller's algorithm for Tate/Weil pairing")
    print("  • Final exponentiation in multiplicative group")
    print("  • BSGS in F_p^k* (can use C++ optimizations)")
    
    print(f"\n✗ Complete MOV implementation requires advanced pairing library.")
    
    return None


def main():
    """Main entry point for optimized MOV attack."""
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
    
    print("="*70)
    print("MOV Attack - Optimized with C++")
    print("="*70)
    print(f"Curve: y² = x³ + {a}x + {b} (mod {p})")
    print(f"Base point G: ({G[0]}, {G[1]})") # type: ignore
    print(f"Target point Q: ({Q[0]}, {Q[1]})") # type: ignore
    print(f"Order n: {n}")
    print()
    
    start_time = time.perf_counter()
    d = mov_attack_optimized(curve, G, Q, n, max_k=12)
    elapsed = time.perf_counter() - start_time
    
    print()
    print("="*70)
    if d is not None:
        Q_verify = fast_scalar_mult(d, G, curve)
        verified = (Q_verify == Q)
        
        print(f"Solution: d = {d}")
        print(f"Time: {elapsed:.6f} seconds")
        print(f"Verification: {'PASSED' if verified else 'FAILED'}")
        print("="*70)
        
        if not verified:
            sys.exit(1)
    else:
        print(f"MOV attack demonstration complete")
        print(f"Time: {elapsed:.6f} seconds")
        print()
        print("Performance Note:")
        print(f"  {'C++' if USE_CPP else 'Python'} backend used for scalar operations")
        print("  Pairing computation would benefit from specialized libraries")
        print("  (e.g., RELIC, PBC, or custom optimized implementations)")
        print("="*70)
        sys.exit(1)


if __name__ == "__main__":
    main()
