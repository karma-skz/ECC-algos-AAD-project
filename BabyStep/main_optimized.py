"""Baby-Step Giant-Step ECDLP - Optimized with C++"""
import sys, time, ctypes, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import EllipticCurve, Point, load_input

# Load C++ library
USE_CPP = False
try:
    lib_path = Path(__file__).parent / 'ecc_fast.so'
    if lib_path.exists():
        ecc_lib = ctypes.CDLL(str(lib_path))
        ecc_lib.scalar_mult.argtypes = [ctypes.c_int64]*6 + [ctypes.POINTER(ctypes.c_int64)]*2
        ecc_lib.scalar_mult.restype = ctypes.c_int
        USE_CPP = True
        print("✓ Using C++ optimization")
except: pass

def fast_scalar_mult(k, G, curve):
    if not USE_CPP or G is None: return curve.scalar_multiply(k, G)
    rx, ry = ctypes.c_int64(), ctypes.c_int64()
    valid = ecc_lib.scalar_mult(ctypes.c_int64(k), ctypes.c_int64(G[0]), ctypes.c_int64(G[1]),
                                  ctypes.c_int64(curve.a), ctypes.c_int64(curve.b), ctypes.c_int64(curve.p),
                                  ctypes.byref(rx), ctypes.byref(ry))
    return None if valid == 0 else (rx.value, ry.value)

if __name__ == "__main__":
    p, a, b, G, n, Q = load_input(Path(sys.argv[1]) if len(sys.argv) > 1 else Path('input/test_1.txt'))
    curve = EllipticCurve(a, b, p)
    print(f"Baby-Step Giant-Step ECDLP: p={p}, n={n}")
    
    start = time.time()
    m = int(math.ceil(math.sqrt(n)))
    
    # Baby steps: compute jG for j=0..m-1
    print(f"Building table (m={m})...")
    baby_table = {}
    R = None
    for j in range(m):
        if j % 1000 == 0: print(f"  Baby step {j}/{m}...", end='\r', flush=True)
        baby_table[R] = j
        R = curve.add(R, G) if R else G
    
    # Giant steps
    print(f"\nSearching...")
    mG = fast_scalar_mult(m, G, curve)
    neg_mG = (mG[0], (-mG[1]) % p) if mG else None
    gamma = Q
    
    for i in range(m):
        if i % 1000 == 0: print(f"  Giant step {i}/{m}...", end='\r', flush=True)
        if gamma in baby_table:
            j = baby_table[gamma]
            d = (i * m + j) % n
            if fast_scalar_mult(d, G, curve) == Q:
                elapsed = time.time() - start
                print(f"\n✓ Solution: d = {d}")
                print(f"Time: {elapsed:.6f}s")
                print("Verification: PASSED")
                sys.exit(0)
        gamma = curve.add(gamma, neg_mG)
    
    print("\n✗ No solution found")
    sys.exit(1)
