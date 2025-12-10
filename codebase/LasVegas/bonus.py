"""
Las Vegas with Approximate Hint - BONUS
Adaptation: Gaussian sampling around approximate key.
"""
import sys, time, random, ctypes, argparse
from pathlib import Path

# --- BOILERPLATE: PATHS & C++ ---
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
    valid = ecc_lib.scalar_mult(k % curve.p, G[0], G[1], curve.a, curve.b, curve.p, ctypes.byref(rx), ctypes.byref(ry)) #type:ignore
    return (rx.value, ry.value) if valid else None

def solve_gaussian(curve, G, Q, n, hint, sigma):
    start = time.perf_counter()
    # INCREASED LIMIT: 50k -> 500k to catch 3-sigma tail events
    limit = 500000 
    
    for i in range(limit):
        # Sample from Gaussian distribution
        offset = int(random.gauss(0, sigma))
        cand = (hint + offset) % n
        
        if fast_mult(cand, G, curve) == Q:
            return cand, i+1, time.perf_counter() - start
            
    return None, limit, time.perf_counter() - start

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Path to test case file")
    parser.add_argument("--approx-error", type=int, help="Approximate error margin (sigma*3)")
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
        else: d_real = 12345
    except: d_real = 12345

    # --- MODE 1: Specific Error Test ---
    if args.approx_error is not None:
        error = args.approx_error
        # Hint is exactly 'error' away to test worst case, or random
        hint = d_real + random.randint(-error, error)
        sigma = max(1, error / 3.0)
        
        d, tries, t = solve_gaussian(curve, G, Q, n, hint, sigma)
        
        print_bonus_result("LasVegas", "success" if d == d_real else "fail", t, tries, {"approx_error": error})
        return

    # --- MODE 2: Default Demo ---
    print(f"\n{'='*60}")
    print(f"LAS VEGAS: GAUSSIAN SAMPLING")
    print(f"{'='*60}")
    
    # 1000 error margin (Scenario: Power Analysis leak)
    error = 1000
    # Simulate a hint that is 'error' away (worst case testing) or random
    hint = d_real + random.randint(-error, error)
    
    # Sigma = error / 3 covers 99.7% of cases
    sigma = error / 3.0
    
    print(f"Hint: {hint} (Error +/- {error})")
    print(f"Sigma: {sigma:.1f}")
    
    d, tries, t = solve_gaussian(curve, G, Q, n, hint, sigma)
    
    if d == d_real:
        print(f"Found: {d} in {tries} tries ({t:.6f}s)")
        print(f"Result: [SUCCESS]")
    else:
        print(f"Result: [FAILED] (Time: {t:.6f}s, Tries: {tries})")
        print("Note: Failed to find key in tail of distribution.")

    print_bonus_result("LasVegas", "success" if d == d_real else "fail", t, tries, {"sigma": sigma})

if __name__ == "__main__":
    main()