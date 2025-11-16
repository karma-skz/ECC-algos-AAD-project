#!/usr/bin/env python3
"""
Bonus ECDLP Algorithm Comparison Script

Tests bonus implementations (with key leakage) on test cases.
- 10-40 bits: All 3 algorithms (BruteForce, BabyStep, PohligHellman)
- 40-50 bits: Only BruteForce, BabyStep, PohligHellman (if feasible)

PollardRho and LasVegas bonus not included as they're infeasible beyond 18 bits.
"""

import sys
import subprocess
import time
from pathlib import Path

# Only deterministic algorithms for bonus (probabilistic ones not feasible)
BONUS_ALGORITHMS = ['BruteForce', 'BabyStep', 'PohligHellman']

def discover_case_files(bit_length: int):
    """Return sorted list of case files for given bitsize."""
    cases_dir = Path('test_cases') / f'{bit_length:02d}bit'
    if not cases_dir.exists():
        return []
    import re
    def case_key(p: Path):
        m = re.search(r'case_(\d+)\.txt$', p.name)
        return int(m.group(1)) if m else 0
    return sorted(cases_dir.glob('case_*.txt'), key=case_key)

def test_bonus_algorithm(algo, bit_length):
    """Test bonus implementation on all available test cases for a bit length."""
    results = []
    case_files = discover_case_files(bit_length)
    for test_file in case_files:
        
        script = Path(algo) / 'bonus.py'
        if not script.exists():
            results.append((None, None, None))
            continue
        
        try:
            start = time.time()
            result = subprocess.run(
                ['python3', str(script), str(test_file)],
                capture_output=True,
                text=True,
                timeout=30  # 30 seconds timeout
            )
            elapsed = time.time() - start
            
            # Extract speedup and leak info
            speedup = None
            leak_bits = None
            
            if result.returncode == 0 and ("PASSED" in result.stdout or "✓ Found:" in result.stdout):
                # Try to extract speedup
                import re
                speedup_match = re.search(r'Speedup.*?(\d+\.?\d*)x', result.stdout)
                if speedup_match:
                    speedup = float(speedup_match.group(1))
                
                # Extract leaked bits
                leak_match = re.search(r'(\d+) bits? leaked', result.stdout, re.IGNORECASE)
                if leak_match:
                    leak_bits = int(leak_match.group(1))
                
                results.append((elapsed, speedup, leak_bits))
            elif "requires exact point order" in result.stdout:
                results.append((-1, None, None))  # N/A marker
            else:
                results.append((None, None, None))
        except subprocess.TimeoutExpired:
            results.append((None, None, None))
        except Exception as e:
            results.append((None, None, None))
    
    return results

def format_bonus_results(results):
    """Format bonus test results showing all cases and average."""
    valid_results = [(t, s) for t, s, l in results if t is not None and t > 0]
    
    # Format individual times
    time_parts = []
    for time_val, speedup_val, leak_val in results:
        if time_val and time_val > 0:
            time_str = f"{time_val*1000:5.1f}ms"
        elif time_val is None:
            time_str = "FAIL "
        else:
            time_str = " N/A "
        time_parts.append(time_str)
    
    time_str = " | ".join(time_parts)
    
    # Calculate average
    if valid_results:
        times = [t for t, s in valid_results]
        speedups = [s for t, s in valid_results if s is not None]
        
        avg_time = sum(times) / len(times)
        if avg_time < 1:
            avg_str = f"avg={avg_time*1000:.1f}ms"
        else:
            avg_str = f"avg={avg_time:.3f}s"
        
        if speedups:
            avg_speedup = sum(speedups) / len(speedups)
            avg_str += f" [{avg_speedup:.1f}x speedup]"
        
        denom = len(results) if results else 0
        passed = f"({len(valid_results)}/{denom})"
    else:
        avg_str = "FAILED"
        denom = len(results) if results else 0
        passed = f"(0/{denom})"
    
    return f"{time_str} | {avg_str:30s} {passed}"

def main():
    if len(sys.argv) > 1:
        try:
            bit_start = int(sys.argv[1])
            bit_end = int(sys.argv[2]) if len(sys.argv) > 2 else bit_start
        except:
            print("Usage: python3 run_bonus_comparison.py [bit_start] [bit_end]")
            sys.exit(1)
    else:
        bit_start, bit_end = 10, 40
    
    print("=" * 130)
    print(f"ECC ECDLP BONUS Performance Comparison ({bit_start}-{bit_end} bits)")
    print("With Key Leakage Optimization")
    print("=" * 130)
    print("Format: Per-case times | Average (with speedup)")
    print("Note: PohligHellman shows N/A when order factorization fails")
    print("      PollardRho/LasVegas not included (infeasible beyond 18 bits)")
    print("=" * 130)
    
    for bits in range(bit_start, bit_end + 1):
        print(f"\n{bits}-bit:")
        for algo in BONUS_ALGORITHMS:
            results = test_bonus_algorithm(algo, bits)
            formatted = format_bonus_results(results)
            print(f"  {algo:15s}: {formatted}")

if __name__ == "__main__":
    main()
