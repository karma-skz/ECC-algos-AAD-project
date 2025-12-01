"""
Pollard's Kangaroo (Lambda) - BONUS
Adaptation: Solves ECDLP in bounded interval [a, b]. O(sqrt(width)).
"""
import sys, time, math, random, ctypes
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import EllipticCurve, Point, load_input
from utils.bonus_utils import print_bonus_result

USE_CPP = False
ecc_lib = None
lib_paths = [
    Path(__file__).parent.parent / "utils" / "cpp" / "ecc_fast.so",
    Path("ecc_fast.so")
]
for p in lib_paths:
    if p.exists():
        try:
            ecc_lib = ctypes.CDLL(str(p))
            ecc_lib.scalar_mult.argtypes = [ctypes.c_longlong]*6 + [ctypes.POINTER(ctypes.c_longlong)]*2
            ecc_lib.scalar_mult.restype = ctypes.c_int
            ecc_lib.point_add.argtypes = [ctypes.c_longlong]*6 + [ctypes.POINTER(ctypes.c_longlong)]*2
            ecc_lib.point_add.restype = ctypes.c_int
            USE_CPP = True
            break
        except: pass

def fast_mult(k, G, curve):
    if not USE_CPP: return curve.scalar_multiply(k, G)
    rx, ry = ctypes.c_longlong(), ctypes.c_longlong()
    # FIX: Pass k directly, do NOT mod p
    valid = ecc_lib.scalar_mult(k, G[0], G[1], curve.a, curve.b, curve.p, ctypes.byref(rx), ctypes.byref(ry))
    return (rx.value, ry.value) if valid else None

def fast_add(P, Q, curve):
    if not USE_CPP: return curve.add(P, Q)
    rx, ry = ctypes.c_longlong(), ctypes.c_longlong()
    valid = ecc_lib.point_add(P[0], P[1], Q[0], Q[1], curve.a, curve.p, ctypes.byref(rx), ctypes.byref(ry))
    return (rx.value, ry.value) if valid else None

def kangaroo(curve, G, Q, lower, upper):
    width = upper - lower
    m = int(math.sqrt(width))
    k = 32 
    
    # Try multiple attempts with different jump sets to ensure success
    for attempt in range(10):
        # Generate random jumps
        # Use attempt in seed to vary jumps on retry
        random.seed(curve.a ^ curve.b ^ attempt)
        
        jumps = []
        for _ in range(k):
            # Mean step size approx sqrt(width)/2
            val = random.randint(1, int(width**0.5) + 1)
            jumps.append(val)
            
        jump_points = [fast_mult(j, G, curve) for j in jumps]
        
        def get_index(P):
            return (P[0] ^ P[1]) % k
        
        # 1. Tame Kangaroo (Trap)
        tame_pos = fast_mult(upper, G, curve)
        tame_dist = 0
        
        # Walk further to increase trap probability (2.0 * m)
        walk_len = int(2.0 * m)
        for _ in range(walk_len):
            idx = get_index(tame_pos)
            tame_pos = fast_add(tame_pos, jump_points[idx], curve)
            tame_dist += jumps[idx]
            
        trap_pos = tame_pos
        total_trap_scalar = upper + tame_dist
        
        # 2. Wild Kangaroo
        wild_pos = Q
        wild_dist = 0
        
        # Limit: width + tame_dist + margin
        limit = width + tame_dist + 2000
        
        while wild_dist < limit:
            if wild_pos == trap_pos:
                return total_trap_scalar - wild_dist
                
            idx = get_index(wild_pos)
            wild_pos = fast_add(wild_pos, jump_points[idx], curve)
            wild_dist += jumps[idx]
            
    return None

def main():
    if len(sys.argv) < 2: return
    p, a, b, G, n, Q = load_input(Path(sys.argv[1]))
    curve = EllipticCurve(a, b, p)
    try:
        num = Path(sys.argv[1]).stem.split('_')[1]
        with open(Path(sys.argv[1]).parent / f"answer_{num}.txt") as f: d_real = int(f.read())
    except: d_real = n // 2

    print(f"\n{'='*70}")
    print(f"POLLARD'S KANGAROO (LAMBDA) DEMO")
    print(f"{'='*70}")
    
    # Smart Interval Selection
    # If curve is tiny (20 bit), 100k width is larger than n.
    # We clamp width to be meaningful (e.g. 1/4 of n)
    target_width = 100000
    if target_width > n:
        target_width = n // 2
        
    lower = max(1, d_real - target_width//2)
    upper = lower + target_width
    
    print(f"Interval: [{lower}, {upper}] (Width: {target_width:,})")
    print(f"Secret:   {d_real}")
    
    t0 = time.perf_counter()
    d = kangaroo(curve, G, Q, lower, upper)
    t = time.perf_counter() - t0
    
    print(f"{'-'*70}")
    if d == d_real:
        print(f"Result: [SUCCESS] (d={d})")
    else:
        print(f"Result: [FAILED] (Found {d}, Expected {d_real})")
    print(f"Time:   {t:.6f}s")

    print_bonus_result("PollardRho", "success" if d == d_real else "fail", t, 0, {"interval_width": target_width})

if __name__ == "__main__":
    main()