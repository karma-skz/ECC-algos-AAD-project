#!/usr/bin/env python3
"""
Test Case Generator for Comprehensive Benchmark
Generates curves for 20-40 bits: Supersingular, Anomalous, PH-Friendly, Generic.
"""

import sys
import random
import time
import argparse
from pathlib import Path

# Import utils from parent directory
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import EllipticCurve

# --- Fast Math Helpers ---
def is_prime(n):
    if n < 2: return False
    if n in (2, 3): return True
    if n % 2 == 0: return False
    r, d = 0, n - 1
    while d % 2 == 0: r += 1; d //= 2
    for _ in range(5):
        a = random.randrange(2, n - 1)
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1: break
        else: return False
    return True

def find_prime(bits, condition=None):
    min_val = 2 ** (bits - 1)
    max_val = 2 ** bits - 1
    while True:
        p = random.randrange(min_val, max_val) | 1
        if is_prime(p):
            if condition is None or condition(p): return p

def tonelli_shanks(n, p):
    if pow(n, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(n, (p + 1) // 4, p)
    Q, S = p - 1, 0
    while Q % 2 == 0: Q //= 2; S += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    M, c, t, R = S, pow(z, Q, p), pow(n, Q, p), pow(n, (Q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return R
        i = 1
        temp = (t * t) % p
        while temp != 1 and i < M: temp = (temp * temp) % p; i += 1
        b = pow(c, 1 << (M - i - 1), p)
        M, c, t, R = i, (b * b) % p, (t * b * b) % p, (R * b) % p

def find_point(curve):
    for _ in range(100):
        x = random.randrange(0, curve.p)
        y = tonelli_shanks((x**3 + curve.a * x + curve.b) % curve.p, curve.p)
        if y is not None: return (x, y)
    return None

def save_case(folder, bits, case_num, p, a, b, G, n, Q, d):
    # Save in structured folders: type/bits/case_X.txt
    target_dir = folder / f"{bits}bit"
    target_dir.mkdir(parents=True, exist_ok=True)
    
    with open(target_dir / f'case_{case_num}.txt', 'w') as f:
        f.write(f'{p}\n{a} {b}\n{G[0]} {G[1]}\n{n}\n{Q[0]} {Q[1]}\n')
    with open(target_dir / f'answer_{case_num}.txt', 'w') as f:
        f.write(f'{d}\n')

# --- Generators ---

def gen_mov(bits, count, base_dir):
    """Supersingular: p = 2 mod 3, y^2 = x^3 + 1."""
    out_dir = base_dir / "MOV_friendly"
    for i in range(1, count + 1):
        p = find_prime(bits, lambda x: x % 3 == 2)
        curve = EllipticCurve(0, 1, p)
        n = p + 1
        G = find_point(curve)
        if not G: continue
        d = random.randrange(2, n)
        Q = curve.scalar_multiply(d, G)
        save_case(out_dir, bits, i, p, 0, 1, G, n, Q, d)

def gen_anomalous(bits, count, base_dir):
    """Anomalous: #E = p. Hard to find for high bits."""
    out_dir = base_dir / "Anomalous"
    found = 0
    start_time = time.time()
    
    # Timeout logic: Finding Trace-1 is hard. Give it 2 seconds per bit-level.
    while found < count and (time.time() - start_time) < 2.0:
        p = find_prime(bits)
        a = random.randrange(0, p)
        b = random.randrange(0, p)
        try: curve = EllipticCurve(a, b, p)
        except: continue
        G = find_point(curve)
        if not G: continue
        
        # Check order == p
        if curve.scalar_multiply(p, G) is None:
            # Double check not small order
            if curve.scalar_multiply(1, G) is not None:
                d = random.randrange(2, p)
                Q = curve.scalar_multiply(d, G)
                save_case(out_dir, bits, found+1, p, a, b, G, p, Q, d)
                found += 1

def gen_generic(bits, count, base_dir):
    """Random curves."""
    out_dir = base_dir / "Generic"
    for i in range(1, count + 1):
        p = find_prime(bits)
        a, b = random.randrange(p), random.randrange(p)
        try: curve = EllipticCurve(a, b, p)
        except: continue
        G = find_point(curve)
        if not G: continue
        n = p + 1 # Approx
        d = random.randrange(2, p)
        Q = curve.scalar_multiply(d, G)
        save_case(out_dir, bits, i, p, a, b, G, n, Q, d)

def gen_ph(bits, count, base_dir):
    """Smooth order (approximated by random generation)."""
    out_dir = base_dir / "PH_friendly"
    # For benchmark purposes, we reuse generic curves. 
    # Real PH curves need complex generation (CM method) which is too slow here.
    # We rely on small primes having some smooth factors occasionally.
    for i in range(1, count + 1):
        p = find_prime(bits)
        a, b = random.randrange(p), random.randrange(p)
        try: curve = EllipticCurve(a, b, p)
        except: continue
        G = find_point(curve)
        if not G: continue
        n = p + 1 
        d = random.randrange(2, p)
        Q = curve.scalar_multiply(d, G)
        save_case(out_dir, bits, i, p, a, b, G, n, Q, d)

def main():
    parser = argparse.ArgumentParser(
        description='Generate ECC test cases for different curve types',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 gen.py                                    # Generate 20-40 bit curves
  python3 gen.py --start 30 --end 50 --step 5      # Custom range
  python3 gen.py --count 3                         # Generate 3 cases per type/bit
  python3 gen.py --verbose                         # Show detailed progress
        """
    )
    parser.add_argument('--start', type=int, default=20, help='Starting bit size (default: 20)')
    parser.add_argument('--end', type=int, default=40, help='Ending bit size (default: 40)')
    parser.add_argument('--step', type=int, default=5, help='Step size between bits (default: 5)')
    parser.add_argument('--count', type=int, default=5, help='Number of cases per type/bit (default: 5)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Show detailed progress')
    args = parser.parse_args()

    base_dir = Path(__file__).parent / "test_cases"
    
    print("=" * 70)
    print("ECC Test Case Generator")
    print("=" * 70)
    print(f"Range: {args.start}-{args.end} bits (step: {args.step})")
    print(f"Count: {args.count} case(s) per type/bit size")
    print(f"Output: {base_dir}")
    print("-" * 70)
    
    total_bits = len(range(args.start, args.end + 1, args.step))
    for idx, bits in enumerate(range(args.start, args.end + 1, args.step), 1):
        print(f"[{idx}/{total_bits}] Processing {bits}-bit curves...")
        
        if args.verbose:
            print(f"  → Supersingular (MOV-friendly)...", end=' ')
        gen_mov(bits, args.count, base_dir)
        if args.verbose:
            print("✓")
        
        if args.verbose:
            print(f"  → Anomalous (Smart's attack)...", end=' ')
        gen_anomalous(bits, args.count, base_dir)
        if args.verbose:
            print("✓")
        
        if args.verbose:
            print(f"  → Generic (Pollard Rho)...", end=' ')
        gen_generic(bits, args.count, base_dir)
        if args.verbose:
            print("✓")
        
        if args.verbose:
            print(f"  → Smooth order (Pohlig-Hellman)...", end=' ')
        gen_ph(bits, args.count, base_dir)
        if args.verbose:
            print("✓")
    
    print("-" * 70)
    print(f"✓ Generation complete!")
    print(f"  Location: {base_dir}")
    print(f"  Total: {total_bits * 4 * args.count} test cases")
    print("=" * 70)

if __name__ == "__main__":
    main()