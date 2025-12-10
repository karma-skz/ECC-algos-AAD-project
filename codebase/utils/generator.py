"""
Test case generator for ECDLP problems.
"""

import random
import os
from pathlib import Path
from typing import Optional, Tuple, Dict
from .ecc_utils import EllipticCurve, Point


def is_prime_miller_rabin(n: int, k: int = 12) -> bool:
    """
    Miller-Rabin primality test.
    
    Args:
        n: Number to test for primality
        k: Number of rounds (higher = more accurate)
    
    Returns:
        True if n is probably prime, False if composite
    """
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    
    # Write n-1 as 2^r * d
    d = n - 1
    r = 0
    while d % 2 == 0:
        d >>= 1
        r += 1
    
    # Witness loop
    for _ in range(k):
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


def generate_prime(bits: int = 40, condition: Optional[int] = None) -> int:
    """
    Generate a random prime of specified bit length.
    
    Args:
        bits: Number of bits in the prime
        condition: If specified, generate p such that p % 12 == condition
    
    Returns:
        A prime number
    """
    while True:
        p = random.getrandbits(bits) | 1  # Ensure odd
        
        if condition is not None and p % 12 != condition:
            continue
        
        if p.bit_length() == bits and is_prime_miller_rabin(p):
            return p


def tonelli_shanks(n: int, p: int) -> Optional[int]:
    """
    Tonelli-Shanks algorithm for computing square roots modulo p.
    Optimized for p ≡ 3 (mod 4).
    
    Args:
        n: Number to find square root of
        p: Prime modulus
    
    Returns:
        Square root of n mod p, or None if n is not a quadratic residue
    """
    # Check if n is a quadratic residue
    if pow(n, (p - 1) // 2, p) != 1:
        return None
    
    # Fast path for p ≡ 3 (mod 4)
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
    
    # General Tonelli-Shanks algorithm
    # Find Q and S such that p - 1 = Q * 2^S
    Q = p - 1
    S = 0
    while Q % 2 == 0:
        Q //= 2
        S += 1
    
    # Find a quadratic non-residue
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    
    M = S
    c = pow(z, Q, p)
    t = pow(n, Q, p)
    R = pow(n, (Q + 1) // 2, p)
    
    while t != 1:
        # Find the least i such that t^(2^i) = 1
        i = 1
        temp = (t * t) % p
        while temp != 1:
            temp = (temp * temp) % p
            i += 1
        
        b = pow(c, 1 << (M - i - 1), p)
        M = i
        c = (b * b) % p
        t = (t * c) % p
        R = (R * b) % p
    
    return R


def trial_factor(n: int) -> Dict[int, int]:
    """
    Factor n using trial division.
    
    Args:
        n: Number to factor
    
    Returns:
        Dictionary mapping prime factors to their exponents
    """
    factors = {}
    d = 2
    
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    
    return factors


def has_full_order(curve: EllipticCurve, G: Point, n: int) -> bool:
    """
    Check if point G has order n on the curve.
    
    Args:
        curve: Elliptic curve
        G: Point to check
        n: Expected order
    
    Returns:
        True if G has order n
    """
    factors = trial_factor(n)
    
    for prime in factors.keys():
        # Check that (n/prime)*G != O
        if curve.scalar_multiply(n // prime, G) is None:
            return False
    
    return True


def generate_test_case(output_dir: Path, case_num: int, secret: Optional[int] = None, bits: int = 40):
    """
    Generate a single ECDLP test case with diverse parameters.
    
    Args:
        output_dir: Directory to save test case files
        case_num: Test case number
        secret: Secret discrete log value (if None, randomly chosen)
        bits: Bit length of prime p
    """
    # Generate prime p ≡ 11 (mod 12) for fast Tonelli
    p = generate_prime(bits=bits, condition=11)
    
    # Generate diverse curve parameters
    # Avoid singular curves: 4a³ + 27b² ≠ 0 (mod p)
    max_curve_attempts = 100
    for _ in range(max_curve_attempts):
        # Generate random a, b with different strategies per case
        if case_num == 1:
            a, b = 0, 1  # Simple curve
        elif case_num == 2:
            a, b = random.randrange(1, 10), random.randrange(1, 10)  # Small coefficients
        elif case_num == 3:
            a, b = random.randrange(p // 2, p), random.randrange(p // 2, p)  # Large coefficients
        elif case_num == 4:
            a = random.randrange(0, p)
            b = random.randrange(0, p)  # Random coefficients
        else:
            a = random.randrange(0, p)
            b = random.randrange(0, p)
        
        # Check non-singularity: 4a³ + 27b² ≠ 0 (mod p)
        discriminant = (4 * pow(a, 3, p) + 27 * pow(b, 2, p)) % p
        if discriminant == 0:
            continue
        
        try:
            curve = EllipticCurve(a, b, p)
            break
        except ValueError:
            continue
    else:
        # Couldn't find valid curve, try again
        return generate_test_case(output_dir, case_num, secret, bits)
    
    # Estimate order using Hasse's theorem: |n - (p+1)| <= 2√p
    # For simplicity, use n ≈ p + 1 and find actual order
    n = p + 1
    
    # Find a point of full order
    max_attempts = 1000
    for _ in range(max_attempts):
        x = random.randrange(0, p)
        rhs = (pow(x, 3, p) + a * x + b) % p
        y = tonelli_shanks(rhs, p)
        
        if y is None:
            continue
        
        G = (x, y)
        if curve.is_on_curve(G) and has_full_order(curve, G, n):
            break
    else:
        # Couldn't find full-order point, try new curve
        return generate_test_case(output_dir, case_num, secret, bits)
    
    # Choose secret value strategically for different cases
    if secret is None:
        sqrt_n = int(n ** 0.5)
        if case_num == 1:
            # Small d < sqrt(n) / 2 - good for brute force
            d = random.randrange(max(1, sqrt_n // 4), sqrt_n // 2)
        elif case_num == 2:
            # Medium d around sqrt(n) / 2 - still manageable for brute force
            d = random.randrange(sqrt_n // 2, (3 * sqrt_n) // 4)
        elif case_num == 3:
            # d around sqrt(n) - competitive between brute force and BSGS
            d = random.randrange((3 * sqrt_n) // 4, sqrt_n + sqrt_n // 10)
        elif case_num == 4:
            # d slightly above sqrt(n) - BSGS should win
            d = random.randrange(sqrt_n + sqrt_n // 10, sqrt_n + sqrt_n // 4)
        else:
            # d well above sqrt(n) - BSGS clear winner
            d = random.randrange(sqrt_n + sqrt_n // 4, sqrt_n + sqrt_n // 2)
    else:
        d = secret % n
    
    # Compute Q = d*G
    Q = curve.scalar_multiply(d, G)
    
    if Q is None:
        # Shouldn't happen, but regenerate if it does
        return generate_test_case(output_dir, case_num, secret, bits)
    
    # Write test case
    testcase_path = output_dir / f"testcase_{case_num}.txt"
    with testcase_path.open('w') as f:
        f.write(f"{p}\n")
        f.write(f"{a} {b}\n")
        f.write(f"{G[0]} {G[1]}\n")
        f.write(f"{n}\n")
        f.write(f"{Q[0]} {Q[1]}\n")
    
    # Write answer
    answer_path = output_dir / f"answer_{case_num}.txt"
    with answer_path.open('w') as f:
        f.write(f"{d}\n")
    
    sqrt_n = int(n ** 0.5)
    d_pos = "small" if d < sqrt_n else "large"
    print(f"Generated test case {case_num}: p={p}, n={n}, d={d} ({d_pos}, √n≈{sqrt_n})")
    print(f"  Curve: y² = x³ + {a}x + {b} (mod {p})")


def generate_test_suite(output_dir: Path, num_cases: int = 5, bits: int = 40):
    """
    Generate a suite of ECDLP test cases with diverse parameters.
    
    Args:
        output_dir: Directory to save test cases
        num_cases: Number of test cases to generate
        bits: Bit length of primes
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Generating {num_cases} diverse test cases...")
    print(f"Prime bit length: {bits}")
    print("="*60)
    
    for i in range(1, num_cases + 1):
        generate_test_case(output_dir, i, secret=None, bits=bits)
        print()
    
    print("="*60)
    print(f"Test suite generated successfully in {output_dir}")


if __name__ == "__main__":
    # Generate test suite in the input directory
    script_dir = Path(__file__).parent.parent
    input_dir = script_dir / "input"
    
    generate_test_suite(input_dir, num_cases=5, bits=40)
