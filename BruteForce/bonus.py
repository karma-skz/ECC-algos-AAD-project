"""
Brute Force with Key Leakage - BONUS
Adaptation: Constrains search to specific bits or intervals.
"""
import sys
import time
import ctypes
import argparse
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
            USE_CPP = True
            break
        except: pass

def fast_mult(k, G, curve):
    if not USE_CPP: return curve.scalar_multiply(k, G)
    rx, ry = ctypes.c_longlong(), ctypes.c_longlong()
    valid = ecc_lib.scalar_mult(k % curve.p, G[0], G[1], curve.a, curve.b, curve.p, ctypes.byref(rx), ctypes.byref(ry))
    return (rx.value, ry.value) if valid else None

def solve_lsb(curve, G, Q, n, leak, bits):
    """Search d = leak (mod 2^bits). Step size = 2^bits."""
    step = 1 << bits
    curr_d = leak
    curr_P = fast_mult(curr_d, G, curve)
    
    # Pre-calc step point
    StepG = fast_mult(step, G, curve)
    
    start = time.perf_counter()
    while curr_d < n:
        if curr_P == Q:
            return curr_d
        curr_P = curve.add(curr_P, StepG)
        curr_d += step
        # Timeout after 5s
        if time.perf_counter() - start > 5: return None 
    return None

def solve_interval(curve, G, Q, lower, upper):
    """Search d in [lower, upper]."""
    curr_P = fast_mult(lower, G, curve)
    start = time.perf_counter()
    
    for d in range(lower, upper + 1):
        if curr_P == Q:
            return d
        curr_P = curve.add(curr_P, G)
        if time.perf_counter() - start > 5: return None
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
        # Try to find answer file
        case_path = Path(args.file)
        # Handle both test_X.txt and testcase_X.txt patterns
        name_parts = case_path.stem.split('_')
        num = name_parts[-1]
        
        ans_file = case_path.parent / f"answer_{num}.txt"
        if not ans_file.exists():
             ans_file = case_path.parent / "answer.txt"
             
        if ans_file.exists():
            with open(ans_file) as f: d_real = int(f.read())
        else:
            # Fallback for demo if no answer file
            d_real = 12345
    except: d_real = 12345

    # --- MODE 1: Specific Leak Test (Called by Analysis Script) ---
    if args.leak_bits is not None:
        bits_to_leak = args.leak_bits
        leak = d_real & ((1<<bits_to_leak)-1)
        space = n // (1<<bits_to_leak)
        
        t0 = time.perf_counter()
        d = solve_lsb(curve, G, Q, n, leak, bits_to_leak)
        t = time.perf_counter() - t0
        
        print_bonus_result("BruteForce", "success" if d == d_real else "fail", t, space, {"leaked_bits": bits_to_leak})
        return

    # --- MODE 2: Specific Interval Test (Called by Analysis Script) ---
    if args.interval_width is not None:
        width = args.interval_width
        lower = max(1, d_real - width//2)
        t0 = time.perf_counter()
        d = solve_interval(curve, G, Q, lower, lower + width)
        t = time.perf_counter() - t0
        
        print_bonus_result("BruteForce", "success" if d == d_real else "fail", t, width, {"interval_width": width})
        return

    # --- MODE 3: Default Demo (No args) ---
    print(f"\n{'='*70}")
    print(f"BRUTE FORCE ADAPTATIONS (d={d_real})")
    print(f"{'='*70}")
    print(f"{'SCENARIO':<25} | {'SPACE':<12} | {'TIME':<10} | {'SPEEDUP'}")
    print(f"{'-'*70}")

    # Baseline (Mocked for large curves)
    t0 = time.perf_counter()
    if n < 2000000: 
        solve_interval(curve, G, Q, 1, n)
        base_t = max(0.000001, time.perf_counter() - t0)
    else:
        # Theoretical time: n / 1,000,000 ops/sec
        base_t = n / 1000000.0 
        print(f"{'Standard':<25} | {n:<12,} | {base_t:.1f}s (Est)| 1.0x")

    # 1. Interval
    width = 10000
    lower = max(1, d_real - width//2)
    t0 = time.perf_counter()
    d = solve_interval(curve, G, Q, lower, lower + width)
    t = time.perf_counter() - t0
    
    speedup_str = f"{base_t/t:.1e}x" if t > 0 else "Inf"
    print(f"{f'Interval (size {width})':<25} | {width:<12,} | {t:.6f}s   | {speedup_str}")

    # 2. LSB (Dynamic)
    # Leak enough bits so remaining work is ~20 bits (1M ops)
    total_bits = n.bit_length()
    bits_to_leak = max(8, total_bits - 20)
    
    leak = d_real & ((1<<bits_to_leak)-1)
    space = n // (1<<bits_to_leak)
    
    t0 = time.perf_counter()
    d = solve_lsb(curve, G, Q, n, leak, bits_to_leak)
    t = time.perf_counter() - t0
    
    speedup_str = f"{base_t/t:.1e}x" if t > 0 else "Inf"
    print(f"{f'Known {bits_to_leak} LSBs':<25} | {space:<12,} | {t:.6f}s   | {speedup_str}")

    print_bonus_result("BruteForce", "success" if d == d_real else "fail", t, space, {"speedup": speedup_str, "leaked_bits": bits_to_leak})

if __name__ == "__main__":
    main()