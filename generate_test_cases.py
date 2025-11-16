#!/usr/bin/env python3
"""
ECC ECDLP Test Case Generator

Usage:
    - Single bitsize (default behavior):
            python3 generate_test_cases.py <k>
                where k = bit length (10-50)
                generates 5 cases by default (unchanged)

    - Range with optional cases per bitsize:
            python3 generate_test_cases.py <start_bit> <end_bit> [cases_per_bit]
            python3 generate_test_cases.py --start <s> --end <e> [--cases <n>]
            python3 generate_test_cases.py <s> <e> --cases <n>
"""

import sys
import random
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils import EllipticCurve

def is_prime(n):
    """Miller-Rabin primality test."""
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False
    
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2
    
    for _ in range(5):
        a = random.randrange(2, n - 1)
        x = pow(a, d, n)
        
        if x == 1 or x == n - 1:
            continue
        
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def find_random_prime(bits, seed_offset=0):
    """Find a random prime with specified bit length."""
    random.seed(bits * 1000 + seed_offset)
    
    min_val = 2 ** (bits - 1)
    max_val = 2 ** bits - 1
    
    for _ in range(10000):
        candidate = random.randrange(min_val, max_val) | 1
        if is_prime(candidate):
            return candidate
    
    candidate = min_val + 1
    while candidate < max_val:
        if is_prime(candidate):
            return candidate
        candidate += 2
    
    return None

def generate_curve_params(p, case_num):
    """Generate varied curve parameters."""
    random.seed(p + case_num)
    
    # Use y^2 = x^3 + ax + b with small a, b
    # Keep it simple for better compatibility
    curves = [
        (0, 1),   # y^2 = x^3 + 1
        (0, 7),   # y^2 = x^3 + 7 (secp256k1-like)
        (1, 1),   # y^2 = x^3 + x + 1
        (2, 3),   # y^2 = x^3 + 2x + 3
        (1, 0),   # y^2 = x^3 + x
    ]
    
    a, b = curves[(case_num - 1) % len(curves)]
    
    # Ensure valid discriminant
    discriminant = (4 * a**3 + 27 * b**2) % p
    if discriminant == 0:
        a, b = (0, 7)  # Fallback
    
    return a, b

def tonelli_shanks(n, p):
    """Compute square root of n modulo prime p using Tonelli-Shanks algorithm."""
    if pow(n, (p - 1) // 2, p) != 1:
        return None  # n is not a quadratic residue
    
    # Find Q and S such that p - 1 = Q * 2^S
    Q, S = p - 1, 0
    while Q % 2 == 0:
        Q //= 2
        S += 1
    
    # Find a quadratic non-residue z
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    
    M = S
    c = pow(z, Q, p)
    t = pow(n, Q, p)
    R = pow(n, (Q + 1) // 2, p)
    
    while True:
        if t == 0:
            return 0
        if t == 1:
            return R
        
        # Find the least i such that t^(2^i) = 1
        i = 1
        temp = (t * t) % p
        while temp != 1 and i < M:
            temp = (temp * temp) % p
            i += 1
        
        b = pow(c, 1 << (M - i - 1), p)
        M = i
        c = (b * b) % p
        t = (t * c) % p
        R = (R * b) % p

def find_generator_point(curve, p, seed_val=0):
    """Find a valid point on the curve using Tonelli-Shanks for sqrt."""
    random.seed(p + seed_val * 1000)
    
    # Try more attempts for larger primes
    max_tries = min(10000, p)
    
    for _ in range(max_tries):
        x = random.randint(0, p - 1)
        y_sq = (x**3 + curve.a * x + curve.b) % p
        
        # Use Tonelli-Shanks to find square root
        y = tonelli_shanks(y_sq, p)
        if y is not None:
            return (x, y)
    
    return None

def generate_test_cases_for_bits(k, num_cases=5):
    """Generate test cases for k-bit ECDLP problems.

    Defaults to 5 cases per bitsize to preserve original behavior.
    """
    if k < 10 or k > 50:
        print(f"Error: Bit length must be between 10 and 50")
        return 0
    
    print(f"\nGenerating {num_cases} test cases for {k}-bit ECDLP...")
    print("=" * 70)
    
    success_count = 0
    
    for case_num in range(1, num_cases + 1):
        p = find_random_prime(k, seed_offset=case_num)
        if not p:
            print(f"  ✗ Case {case_num}: No prime found")
            continue
        
        a, b = generate_curve_params(p, case_num)
        curve = EllipticCurve(a, b, p)
        
        G = find_generator_point(curve, p, seed_val=case_num)
        if not G:
            print(f"  ✗ Case {case_num}: No generator found")
            continue
        
        # Use small discrete log for solvability
        # For smaller bit sizes, keep secrets small for testing
        # For larger bit sizes (>20), use moderately sized secrets to test algorithms properly
        if k <= 18:
            # Small secrets for quick testing of all 5 algorithms
            d = 50 + case_num * 50 + k * 5
            max_d = min(5000, 2**(k-3))
        elif k <= 25:
            # Medium secrets - test BSGS properly
            d = 1000 + case_num * 500 + k * 100
            max_d = min(50000, 2**(k-5))
        else:
            # Larger secrets for realistic testing (but still solvable)
            d = 5000 + case_num * 2000 + k * 500
            max_d = min(500000, 2**(k-8))
        
        d = d % max_d
        if d < 10:
            d += 10
        
        Q = curve.scalar_multiply(d, G)
        
        # Use p+1 as order approximation
        n = p + 1
        
        test_dir = Path(__file__).parent / 'test_cases' / f'{k:02d}bit'
        test_dir.mkdir(parents=True, exist_ok=True)
        
        filename = test_dir / f'case_{case_num}.txt'
        with open(filename, 'w') as f:
            f.write(f'{p}\n')
            f.write(f'{a} {b}\n')
            f.write(f'{G[0]} {G[1]}\n')
            f.write(f'{n}\n')
            f.write(f'{Q[0]} {Q[1]}\n') # type: ignore
        
        # Also save the answer (private key) for bonus implementations
        answer_file = test_dir / f'answer_{case_num}.txt'
        with open(answer_file, 'w') as f:
            f.write(f'{d}\n')
        
        print(f"  ✓ Case {case_num}: p={p} ({k} bits), a={a}, b={b}, d={d}")
        success_count += 1
    
    print("=" * 70)
    print(f"✓ Generated {success_count}/{num_cases} test cases")
    return success_count

def main():
    args = sys.argv[1:]

    # Minimal argument parsing supporting positional and flags while
    # preserving the original defaults (single k -> 5 cases).
    def print_usage_and_exit():
        print("Usage:")
        print("  python3 generate_test_cases.py <k>")
        print("  python3 generate_test_cases.py <start_bit> <end_bit> [cases_per_bit]")
        print("  python3 generate_test_cases.py --start <s> --end <e> [--cases <n>]")
        print("Constraints: bit length must be between 10 and 50")
        sys.exit(1)

    # Extract optional --start/--end/--cases flags if provided
    start_flag = None
    end_flag = None
    cases_flag = None
    i = 0
    positional = []
    while i < len(args):
        if args[i] in ("--start", "-s"):
            if i + 1 >= len(args):
                print("Error: --start requires a value")
                print_usage_and_exit()
            start_flag = args[i + 1]
            i += 2
        elif args[i] in ("--end", "-e"):
            if i + 1 >= len(args):
                print("Error: --end requires a value")
                print_usage_and_exit()
            end_flag = args[i + 1]
            i += 2
        elif args[i] in ("--cases", "-c"):
            if i + 1 >= len(args):
                print("Error: --cases requires a value")
                print_usage_and_exit()
            cases_flag = args[i + 1]
            i += 2
        else:
            positional.append(args[i])
            i += 1

    # Determine start, end, and cases from provided inputs
    start_bit = None
    end_bit = None
    cases_per_bit = None

    # First, if flags present, use them
    if start_flag is not None or end_flag is not None:
        if start_flag is None or end_flag is None:
            print("Error: --start and --end must be used together")
            print_usage_and_exit()
        try:
            start_bit = int(start_flag) # type: ignore
            end_bit = int(end_flag) # type: ignore
        except ValueError:
            print("Error: --start and --end must be integers")
            sys.exit(1)
    
    # Next, handle positional variants
    if len(positional) == 1 and start_bit is None:
        # Original mode: single k
        try:
            k = int(positional[0])
        except ValueError:
            print("Error: <k> must be an integer")
            sys.exit(1)
        start_bit = k
        end_bit = k
    elif len(positional) == 2 and start_bit is None:
        try:
            start_bit = int(positional[0])
            end_bit = int(positional[1])
        except ValueError:
            print("Error: <start_bit> and <end_bit> must be integers")
            sys.exit(1)
    elif len(positional) == 3 and start_bit is None:
        try:
            start_bit = int(positional[0])
            end_bit = int(positional[1])
            cases_per_bit = int(positional[2])
        except ValueError:
            print("Error: <start_bit> <end_bit> [cases_per_bit] must be integers")
            sys.exit(1)
    elif len(positional) > 3:
        print_usage_and_exit()

    # If neither flags nor positional provided for start/end, show usage
    if start_bit is None or end_bit is None:
        print_usage_and_exit()
    # Type narrowing for linters/type checkers
    assert isinstance(start_bit, int) and isinstance(end_bit, int)

    # If cases specified via flag, it overrides positional third arg
    if cases_flag is not None:
        try:
            cases_per_bit = int(cases_flag)
        except ValueError:
            print("Error: --cases must be an integer")
            sys.exit(1)

    if cases_per_bit is None:
        cases_per_bit = 5  # Preserve original default

    # Validate ranges
    if start_bit < 10 or end_bit > 50 or start_bit > end_bit:
        print("Error: Bit length range must be 10 <= start <= end <= 50")
        sys.exit(1)

    print("=" * 70)
    print(f"ECC ECDLP Test Case Generator")
    print("=" * 70)
    print(f"Bit lengths: {start_bit} to {end_bit} (inclusive)")
    print(f"Cases per bitsize: {cases_per_bit}")

    total_bits = 0
    for k in range(start_bit, end_bit + 1):
        count = generate_test_cases_for_bits(k, num_cases=cases_per_bit)
        if count > 0:
            print(f"  Location: test_cases/{k:02d}bit/")
        total_bits += 1

    # If at least one bitsize produced at least one case, consider success
    # Note: generate_test_cases_for_bits already prints per-bitsize results
    print("\n✓ Done generating test cases across range.")

if __name__ == "__main__":
    main()
