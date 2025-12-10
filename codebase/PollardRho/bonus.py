"""
Pollard's Kangaroo (Lambda) - BONUS
Adaptation: Solves ECDLP in bounded interval [a, b]. O(sqrt(width)).
"""
import sys, time, math, random, ctypes, argparse
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
    valid = ecc_lib.scalar_mult(k, G[0], G[1], curve.a, curve.b, curve.p, ctypes.byref(rx), ctypes.byref(ry)) #type:ignore
    return (rx.value, ry.value) if valid else None

def fast_add(P, Q, curve):
    if not USE_CPP: return curve.add(P, Q)
    rx, ry = ctypes.c_longlong(), ctypes.c_longlong()
    valid = ecc_lib.point_add(P[0], P[1], Q[0], Q[1], curve.a, curve.p, ctypes.byref(rx), ctypes.byref(ry)) #type:ignore
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

def solve_lsb(curve, G, Q, n, leak, bits):
    """
    Solves d = x * 2^bits + leak using Kangaroo on the transformed problem.
    """
    # 1. Compute Q'
    leakG = fast_mult(leak, G, curve)
    Q_prime = curve.add(Q, curve.negate(leakG))
    
    # 2. Compute G'
    factor = 1 << bits
    G_prime = fast_mult(factor, G, curve)
    
    # 3. Solve smaller DLP: Q' = x * G'
    limit = n // factor
    
    # Use Kangaroo on [0, limit]
    # Note: Kangaroo expects (curve, G, Q, lower, upper)
    return kangaroo(curve, G_prime, Q_prime, 0, limit)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Path to test case file")
    parser.add_argument("--interval-width", type=int, help="Width of interval to search")
    parser.add_argument("--leak-bits", type=int, help="Number of LSB bits to leak")
    args = parser.parse_args()

    p, a, b, G, n, Q = load_input(Path(args.file))
    curve = EllipticCurve(a, b, p)
    try:
        case_path = Path(args.file)
        name_parts = case_path.stem.split('_')
        num = name_parts[-1]
        ans_file = case_path.parent / f"answer_{num}.txt"
        if not ans_file.exists(): ans_file = case_path.parent / "answer.txt"
        
        if ans_file.exists():
            with open(ans_file) as f: d_real = int(f.read())
        else:
            d_real = n // 2
    except: d_real = n // 2

    # --- MODE 1: Specific Leak Test ---
    if args.leak_bits is not None:
        bits_to_leak = args.leak_bits
        leak = d_real & ((1<<bits_to_leak)-1)
        
        t0 = time.perf_counter()
        # Result x is the upper part, so d = x * 2^bits + leak
        x = solve_lsb(curve, G, Q, n, leak, bits_to_leak)
        t = time.perf_counter() - t0
        
        d = (x * (1<<bits_to_leak) + leak) if x is not None else None
        
        print_bonus_result("PollardRho", "success" if d == d_real else "fail", t, 0, {"leaked_bits": bits_to_leak})
        return

    # --- MODE 2: Specific Interval Test ---
    if args.interval_width is not None:
        target_width = args.interval_width
        lower = max(1, d_real - target_width//2)
        upper = lower + target_width
        
        t0 = time.perf_counter()
        d = kangaroo(curve, G, Q, lower, upper)
        t = time.perf_counter() - t0
        
        print_bonus_result("PollardRho", "success" if d == d_real else "fail", t, 0, {"interval_width": target_width})
        return

    # --- MODE 2: Default Demo ---
    print(f"\n{'='*70}")
    print(f"POLLARD'S KANGAROO (LAMBDA) DEMO")
    print(f"{'='*70}")
    
    # Smart Interval Selection
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