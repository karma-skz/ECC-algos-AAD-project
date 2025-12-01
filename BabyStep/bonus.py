"""
BSGS with Interval Leakage - BONUS
Adaptation: Reduces RAM/Time by searching only the effective range.
"""
import sys, time, math, ctypes, argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import EllipticCurve, Point, load_input
from utils.bonus_utils import print_bonus_result

USE_CPP = False
ecc_lib = None
lib_paths = [Path(__file__).parent.parent / "utils" / "cpp" / "ecc_fast.so"]
for p in lib_paths:
    if p.exists():
        ecc_lib = ctypes.CDLL(str(p))
        ecc_lib.scalar_mult.argtypes = [ctypes.c_longlong]*6 + [ctypes.POINTER(ctypes.c_longlong)]*2
        ecc_lib.scalar_mult.restype = ctypes.c_int
        USE_CPP = True
        break

def fast_mult(k, G, curve):
    if not USE_CPP: return curve.scalar_multiply(k, G)
    rx, ry = ctypes.c_longlong(), ctypes.c_longlong()
    # FIX: No % p
    valid = ecc_lib.scalar_mult(k, G[0], G[1], curve.a, curve.b, curve.p, ctypes.byref(rx), ctypes.byref(ry))
    return (rx.value, ry.value) if valid else None

def bsgs_range(curve, G, Q, lower, upper):
    Q_prime = curve.add(Q, curve.negate(fast_mult(lower, G, curve)))
    limit = upper - lower
    m = int(math.ceil(math.sqrt(limit)))
    
    table = {}
    P = None 
    for j in range(m):
        table[P] = j
        P = curve.add(P, G) if P else G
        
    mG = fast_mult(m, G, curve)
    neg_mG = curve.negate(mG)
    curr = Q_prime
    
    for i in range(m + 1):
        if curr in table:
            return lower + (i * m + table[curr])
        curr = curve.add(curr, neg_mG)
    return None

def solve_lsb(curve, G, Q, n, leak, bits):
    """
    Solves d = x * 2^bits + leak.
    Q = (x * 2^bits + leak) * G
    Q - leak*G = x * (2^bits * G)
    Let Q' = Q - leak*G, G' = 2^bits * G
    Solve Q' = x * G' for x in [0, n / 2^bits]
    """
    # 1. Compute Q'
    leakG = fast_mult(leak, G, curve)
    Q_prime = curve.add(Q, curve.negate(leakG))
    
    # 2. Compute G'
    factor = 1 << bits
    G_prime = fast_mult(factor, G, curve)
    
    # 3. Solve smaller DLP
    limit = n // factor
    # We can reuse bsgs_range for 0..limit
    # Note: bsgs_range expects (curve, G, Q, lower, upper)
    # Here we pass G_prime as the base point
    
    # We need to adapt bsgs_range to accept a custom base point G_prime
    # Or just inline the logic here
    m = int(math.ceil(math.sqrt(limit)))
    
    table = {}
    P = None 
    for j in range(m):
        table[P] = j
        P = curve.add(P, G_prime) if P else G_prime
        
    mG = fast_mult(m, G_prime, curve)
    neg_mG = curve.negate(mG)
    curr = Q_prime
    
    for i in range(m + 1):
        if curr in table:
            x = i * m + table[curr]
            return x * factor + leak
        curr = curve.add(curr, neg_mG)
    return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Path to test case file")
    parser.add_argument("--leak-bits", type=int, help="Number of LSB bits to leak")
    parser.add_argument("--interval-width", type=int, help="Width of interval to search")
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
        else: d_real = n // 2
    except: d_real = n // 2

    # --- MODE 1: Specific Leak Test ---
    if args.leak_bits is not None:
        bits_to_leak = args.leak_bits
        leak = d_real & ((1<<bits_to_leak)-1)
        
        t0 = time.perf_counter()
        d = solve_lsb(curve, G, Q, n, leak, bits_to_leak)
        t = time.perf_counter() - t0
        
        # Calculate effective table size for reporting
        remaining_space = n // (1<<bits_to_leak)
        m_opt = int(math.sqrt(remaining_space))
        
        print_bonus_result("BabyStep", "success" if d == d_real else "fail", t, m_opt, {"leaked_bits": bits_to_leak})
        return

    # --- MODE 2: Specific Interval Test ---
    if args.interval_width is not None:
        width = args.interval_width
        low = max(1, d_real - width//2)
        t0 = time.perf_counter()
        d = bsgs_range(curve, G, Q, low, low + width)
        t = time.perf_counter() - t0
        m_opt = int(math.sqrt(width))
        
        print_bonus_result("BabyStep", "success" if d == d_real else "fail", t, m_opt, {"interval_width": width})
        return

    # --- MODE 3: Default Demo ---
    print(f"\n{'='*70}")
    print(f"BSGS ADAPTATIONS (d={d_real})")
    print(f"{'='*70}")
    print(f"{'SCENARIO':<25} | {'TABLE SIZE':<12} | {'TIME':<10}")
    print(f"{'-'*70}")

    m_std = int(math.sqrt(n))
    print(f"{'Standard':<25} | {m_std:<12,} | {'(Baseline)':<10}")

    width = n // 100
    low = max(1, d_real - width//2)
    t0 = time.perf_counter()
    d = bsgs_range(curve, G, Q, low, low + width)
    t = time.perf_counter() - t0
    m_opt = int(math.sqrt(width))
    print(f"{'1% Interval Leak':<25} | {m_opt:<12,} | {t:.6f}s")
    print(f"{'-'*70}")
    print(f"RAM Savings: {m_std/max(1, m_opt):.1f}x less memory needed")

    print_bonus_result("BabyStep", "success" if d == d_real else "fail", t, m_opt, {"ram_savings": f"{m_std/max(1, m_opt):.1f}x"})

if __name__ == "__main__":
    main()