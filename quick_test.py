#!/usr/bin/env python3
"""Quick test of all algorithms on 10-14 bit test cases"""
import subprocess, time, json
from pathlib import Path
from collections import defaultdict

ALGORITHMS = {
    'BruteForce': {'script': 'BruteForce/main_optimized.py', 'timeout': 30},
    'BabyStep': {'script': 'BabyStep/main_optimized.py', 'timeout': 30},
    'PohligHellman': {'script': 'PohligHellman/main_optimized.py', 'timeout': 60},
    'PollardRho': {'script': 'PollardRho/main_optimized.py', 'timeout': 60},
    'LasVegas': {'script': 'LasVegas/main_optimized.py', 'timeout': 60}
}

def run_test(script_path, test_file, timeout):
    start = time.time()
    try:
        result = subprocess.run(['python3', script_path, str(test_file)], timeout=timeout,
                              capture_output=True, text=True, cwd=Path(__file__).parent)
        elapsed = time.time() - start
        success = 'PASSED' in result.stdout
        attempts = result.stdout.count('Attempt') if 'Attempt' in result.stdout else 1
        return {'success': success, 'time': elapsed, 'attempts': attempts,
                'status': 'SUCCESS' if success else 'FAILED'}
    except subprocess.TimeoutExpired:
        return {'success': False, 'time': timeout, 'attempts': None, 'status': 'TIMEOUT'}
    except:
        return {'success': False, 'time': 0, 'attempts': None, 'status': 'ERROR'}

results = defaultdict(lambda: defaultdict(list))
project_root = Path(__file__).parent

print("="*80)
print("Quick Test: 10-14 bit cases (all available cases)")
print("="*80)

for bits in range(10, 15):
    print(f"\n{bits}-bit:")
    test_dir = project_root / 'test_cases' / f'{bits:02d}bit'
    if not test_dir.exists():
        print("  (no cases)")
        continue
    import re
    def case_key(p: Path):
        m = re.search(r'case_(\d+)\.txt$', p.name)
        return int(m.group(1)) if m else 0
    case_files = sorted(test_dir.glob('case_*.txt'), key=case_key)
    for test_file in case_files:
        m = re.search(r'case_(\d+)\.txt$', test_file.name)
        case_label = m.group(1) if m else test_file.name
        print(f"  Case {case_label}:", end=' ')
        for algo_name, config in ALGORITHMS.items():
            script_path = project_root / config['script']
            if not script_path.exists(): continue
            result = run_test(script_path, test_file, config['timeout'])
            results[bits][algo_name].append(result)
            symbol = "✓" if result['success'] else "✗"
            print(f"{algo_name[0]}{symbol}", end=' ')
        print()

print("\n" + "="*80)
print("SUMMARY (Success Rate / Avg Time)")
print("="*80)
print(f"{'Bits':<6} {'BruteForce':<15} {'BabyStep':<15} {'PohligHellman':<18} {'PollardRho':<15} {'LasVegas':<15}")
print("-"*80)

for bits in range(10, 15):
    row = f"{bits:<6}"
    for algo_name in ALGORITHMS.keys():
        if algo_name in results[bits]:
            cases = results[bits][algo_name]
            success_count = sum(1 for c in cases if c['success'])
            times = [c['time'] for c in cases if c['success']]
            avg_time = sum(times)/len(times) if times else 0
            denom = len(cases)
            row += f" {success_count}/{denom} {avg_time:>6.3f}s   "
        else:
            row += f" -/-  ------s   "
    print(row)

print("="*80)
