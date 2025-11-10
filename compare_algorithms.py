#!/usr/bin/env python3
"""
Compare all ECDLP algorithms on the same test case.
"""

import sys
import time
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from utils import EllipticCurve, load_input

# Import algorithm functions
sys.path.insert(0, str(Path(__file__).parent / 'BruteForce'))
from BruteForce.main import brute_force_ecdlp as brute_force

sys.path.insert(0, str(Path(__file__).parent / 'BabyStep'))
from BabyStep.main import bsgs_ecdlp

sys.path.insert(0, str(Path(__file__).parent / 'PohligHellman'))
from PohligHellman.main import pohlig_hellman_ecdlp

sys.path.insert(0, str(Path(__file__).parent / 'PollardRho'))
from PollardRho.main import pollard_rho_ecdlp

sys.path.insert(0, str(Path(__file__).parent / 'LasVegas'))
from LasVegas.main import las_vegas_ecdlp


def run_algorithm(name, algo_func, curve, G, Q, n, **kwargs):
    """Run an algorithm and return results."""
    print(f"\n{'='*60}")
    print(f"Algorithm: {name}")
    print(f"{'='*60}")
    
    try:
        start = time.perf_counter()
        result = algo_func(curve, G, Q, n, **kwargs)
        elapsed = time.perf_counter() - start
        
        # Handle different return types
        if isinstance(result, tuple):
            d, extra = result
        else:
            d = result
            extra = None
        
        if d is not None:
            # Verify
            Q_verify = curve.scalar_multiply(d, G)
            verified = (Q_verify == Q)
            
            print(f"✓ Solution found: d = {d}")
            print(f"  Time: {elapsed:.6f} seconds")
            print(f"  Verification: {'PASSED' if verified else 'FAILED'}")
            if extra is not None:
                print(f"  Extra info: {extra}")
            
            return True, elapsed
        else:
            print(f"✗ No solution found")
            print(f"  Time: {elapsed:.6f} seconds")
            return False, elapsed
    
    except Exception as e:
        print(f"✗ Error: {e}")
        return False, 0.0


def main():
    """Compare all algorithms."""
    # Load test case
    input_path = Path(__file__).parent / 'input' / 'testcase_1.txt'
    
    if not input_path.exists():
        # Try BruteForce input
        input_path = Path(__file__).parent / 'BruteForce' / 'input' / 'test_1.txt'
    
    if not input_path.exists():
        print("Error: No test file found")
        sys.exit(1)
    
    print("Loading test case...")
    p, a, b, G, n, Q = load_input(input_path)
    curve = EllipticCurve(a, b, p)
    
    print(f"\nTest Case:")
    print(f"  Curve: y² = x³ + {a}x + {b} (mod {p})")
    print(f"  G = {G}")
    print(f"  Q = {Q}")
    print(f"  n = {n}")
    
    # Run algorithms
    results = {}
    
    # Brute Force
    success, time_bf = run_algorithm("Brute Force", brute_force, curve, G, Q, n)
    results['Brute Force'] = (success, time_bf)
    
    # Baby-Step Giant-Step
    success, time_bsgs = run_algorithm("Baby-Step Giant-Step", bsgs_ecdlp, curve, G, Q, n)
    results['BSGS'] = (success, time_bsgs)
    
    # Pohlig-Hellman
    success, time_ph = run_algorithm("Pohlig-Hellman", pohlig_hellman_ecdlp, curve, G, Q, n)
    results['Pohlig-Hellman'] = (success, time_ph)
    
    # Pollard Rho
    success, time_rho = run_algorithm("Pollard Rho", pollard_rho_ecdlp, curve, G, Q, n, 
                                     max_steps=1_000_000, partition_m=32)
    results['Pollard Rho'] = (success, time_rho)
    
    # Las Vegas
    success, time_lv = run_algorithm("Las Vegas", las_vegas_ecdlp, curve, G, Q, n, 
                                    n_prime=2, max_attempts=50)
    results['Las Vegas'] = (success, time_lv)
    
    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"{'Algorithm':<25} {'Success':<10} {'Time (s)':<15}")
    print(f"{'-'*60}")
    
    for name, (success, elapsed) in results.items():
        status = "✓" if success else "✗"
        print(f"{name:<25} {status:<10} {elapsed:<15.6f}")
    
    # Find fastest
    successful = [(name, t) for name, (s, t) in results.items() if s]
    if successful:
        fastest = min(successful, key=lambda x: x[1])
        print(f"\nFastest successful algorithm: {fastest[0]} ({fastest[1]:.6f}s)")


if __name__ == "__main__":
    main()
