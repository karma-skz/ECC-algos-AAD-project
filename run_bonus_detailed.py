#!/usr/bin/env python3
"""
Detailed Bonus Test Runner with Comprehensive Logging and Visualization

Runs all bonus implementations and logs:
- All scenarios tested (LSB/MSB/intervals/residues)
- Individual scenario results
- Speedup metrics
- Success/failure rates
- Detailed timing information
- Performance graphs
"""

import subprocess
import sys
import time
import json
import re
from pathlib import Path
from datetime import datetime

ALGORITHMS = ['BruteForce', 'BabyStep', 'PohligHellman']

# Try to import matplotlib
try:
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

def run_bonus_test(algo, test_file):
    """Run bonus script and capture full output."""
    script = Path(algo) / 'bonus.py'
    if not script.exists():
        return None, "Script missing"
    
    try:
        start = time.time()
        result = subprocess.run(
            ['python3', str(script), str(test_file)],
            capture_output=True,
            text=True,
            timeout=60
        )
        elapsed = time.time() - start
        
        return {
            'returncode': result.returncode,
            'stdout': result.stdout,
            'stderr': result.stderr,
            'elapsed': elapsed,
            'success': result.returncode == 0
        }, None
    
    except subprocess.TimeoutExpired:
        return None, "Timeout"
    except Exception as e:
        return None, str(e)

def parse_bonus_output(output):
    """Parse bonus script output to extract scenario details."""
    scenarios = []
    current_scenario = None
    
    lines = output.split('\n')
    for i, line in enumerate(lines):
        # Detect scenario headers
        if 'SCENARIO' in line and ':' in line:
            if current_scenario:
                scenarios.append(current_scenario)
            current_scenario = {
                'name': line.split(':', 1)[1].strip() if ':' in line else line,
                'tests': []
            }
        
        # Detect test results
        elif current_scenario and ('✓ Found:' in line or '✗ Failed' in line):
            test_result = {
                'success': '✓' in line,
                'line': line.strip()
            }
            
            # Look for timing in next few lines
            for j in range(i, min(i+3, len(lines))):
                if 'Time:' in lines[j]:
                    test_result['time'] = lines[j].strip()
                elif 'Reduction:' in lines[j]:
                    test_result['reduction'] = lines[j].strip()
                elif 'Search space:' in lines[j]:
                    test_result['search_space'] = lines[j].strip()
            
            current_scenario['tests'].append(test_result)
        
        # Detect speedup comparisons
        elif 'Speedup:' in line:
            if current_scenario:
                current_scenario['speedup'] = line.strip()
    
    if current_scenario:
        scenarios.append(current_scenario)
    
    return scenarios

def create_log_entry(algo, bit_length, case_num, result, error):
    """Create structured log entry."""
    entry = {
        'algorithm': algo,
        'bit_length': bit_length,
        'case_number': case_num,
        'timestamp': datetime.now().isoformat(),
        'success': False,
        'error': error
    }
    
    if result:
        entry['success'] = result['success']
        entry['elapsed_total'] = result['elapsed']
        entry['scenarios'] = parse_bonus_output(result['stdout'])
        entry['returncode'] = result['returncode']
        
        if result['stderr']:
            entry['stderr'] = result['stderr']
    
    return entry

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 run_bonus_detailed.py <bit_length1> [bit_length2] ...")
        print("Example: python3 run_bonus_detailed.py 10 12 14")
        sys.exit(1)
    
    bit_lengths = [int(x) for x in sys.argv[1:]]
    
    # Create logs directory
    log_dir = Path('bonus_logs')
    log_dir.mkdir(exist_ok=True)
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = log_dir / f'bonus_test_{timestamp}.log'
    json_file = log_dir / f'bonus_test_{timestamp}.json'
    
    all_results = []
    
    print("="*100)
    print(f"DETAILED BONUS IMPLEMENTATION TEST")
    print(f"Logging to: {log_file}")
    print("="*100)
    
    with open(log_file, 'w') as f:
        f.write("="*100 + "\n")
        f.write(f"ECC ECDLP Bonus Implementation Detailed Test Log\n")
        f.write(f"Started: {datetime.now()}\n")
        f.write("="*100 + "\n\n")
        
        for bits in bit_lengths:
            print(f"\n{'='*100}")
            print(f"Testing {bits}-bit cases")
            print('='*100)
            
            f.write(f"\n{'='*100}\n")
            f.write(f"{bits}-BIT TEST CASES\n")
            f.write('='*100 + "\n\n")
            
            # Discover available cases
            cases_dir = Path('test_cases') / f'{bits:02d}bit'
            if not cases_dir.exists():
                continue
            def case_key(p: Path):
                m = re.search(r'case_(\d+)\.txt$', p.name)
                return int(m.group(1)) if m else 0
            case_files = sorted(cases_dir.glob('case_*.txt'), key=case_key)

            for test_file in case_files:
                m = re.search(r'case_(\d+)\.txt$', test_file.name)
                case_label = m.group(1) if m else test_file.name
                
                print(f"\n  Case {case_label}:")
                f.write(f"\n--- Case {case_label} ---\n")
                
                for algo in ALGORITHMS:
                    print(f"    {algo:15s}: ", end='', flush=True)
                    
                    result, error = run_bonus_test(algo, test_file)
                    try:
                        case_num_value = int(case_label)
                    except Exception:
                        case_num_value = 0
                    entry = create_log_entry(algo, bits, case_num_value, result, error)
                    all_results.append(entry)
                    
                    if result and result['success']:
                        scenarios = entry['scenarios']
                        print(f"✓ {len(scenarios)} scenarios, {result['elapsed']:.2f}s total")
                        
                        f.write(f"\n{algo}:\n")
                        f.write(f"  Status: SUCCESS\n")
                        f.write(f"  Total time: {result['elapsed']:.3f}s\n")
                        f.write(f"  Scenarios tested: {len(scenarios)}\n\n")
                        
                        for i, scenario in enumerate(scenarios, 1):
                            f.write(f"  Scenario {i}: {scenario['name']}\n")
                            for j, test in enumerate(scenario['tests'], 1):
                                status = "PASS" if test['success'] else "FAIL"
                                f.write(f"    Test {j}: {status}\n")
                                if 'reduction' in test:
                                    f.write(f"      {test['reduction']}\n")
                                if 'time' in test:
                                    f.write(f"      {test['time']}\n")
                            
                            if 'speedup' in scenario:
                                f.write(f"  {scenario['speedup']}\n")
                            f.write("\n")
                    else:
                        print(f"✗ {error or 'Failed'}")
                        f.write(f"\n{algo}: FAILED - {error or 'Unknown error'}\n")
                
                f.write("\n")
        
        f.write("\n" + "="*100 + "\n")
        f.write(f"Test completed: {datetime.now()}\n")
        f.write("="*100 + "\n")
    
    # Save JSON summary
    with open(json_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    
    print("\n" + "="*100)
    print(f"✓ Test complete!")
    print(f"  Detailed log: {log_file}")
    print(f"  JSON data: {json_file}")
    
    # Summary statistics
    total_tests = len(all_results)
    successful = sum(1 for r in all_results if r['success'])
    
    print(f"\nSummary:")
    print(f"  Total tests: {total_tests}")
    print(f"  Successful: {successful} ({successful*100//total_tests if total_tests > 0 else 0}%)")
    print(f"  Failed: {total_tests - successful}")
    
    # Per-algorithm summary and collect plot data
    print(f"\nPer-algorithm breakdown:")
    plot_data = {algo: {'bits': [], 'times': []} for algo in ALGORITHMS}
    
    for algo in ALGORITHMS:
        algo_results = [r for r in all_results if r['algorithm'] == algo]
        algo_success = sum(1 for r in algo_results if r['success'])
        if algo_results:
            total_scenarios = sum(len(r.get('scenarios', [])) for r in algo_results if r['success'])
            print(f"  {algo:15s}: {algo_success}/{len(algo_results)} successful, {total_scenarios} scenarios tested")
            
            # Collect data for plotting (average time per bit length)
            for bits in bit_lengths:
                bit_results = [r for r in algo_results if r['bit_length'] == bits and r['success']]
                if bit_results:
                    avg_time = sum(r['elapsed_total'] for r in bit_results) / len(bit_results)
                    plot_data[algo]['bits'].append(bits)
                    plot_data[algo]['times'].append(avg_time)
    
    # Generate plots
    if HAS_MATPLOTLIB:
        print("\n" + "="*100)
        print("Generating performance graphs...")
        generate_bonus_plots(plot_data, bit_lengths, log_dir)
        print("✓ Graphs saved!")
    
    print("="*100)

def generate_bonus_plots(plot_data, bit_lengths, log_dir):
    """Generate performance visualization for bonus implementations."""
    
    colors = {
        'BruteForce': '#e74c3c',
        'BabyStep': '#3498db',
        'PohligHellman': '#2ecc71'
    }
    
    # Create plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6)) # type: ignore
    
    # Linear scale
    max_time = 0
    for algo in ALGORITHMS:
        if plot_data[algo]['bits']:
            times = plot_data[algo]['times']
            max_time = max(max_time, max(times))
            ax1.plot(plot_data[algo]['bits'], times, 
                    marker='o', label=f'{algo} (bonus)', color=colors[algo], linewidth=2, markersize=8)
    
    ax1.set_xlabel('Bit Length', fontsize=12)
    ax1.set_ylabel('Average Time (seconds)', fontsize=12)
    ax1.set_title(f'Bonus Implementation Performance - Linear Scale', fontsize=14, fontweight='bold')
    ax1.legend(loc='upper left')
    ax1.grid(True, alpha=0.3)
    
    # Clip for readability
    clip_threshold = max_time * 0.8 if max_time > 5 else None
    if clip_threshold:
        ax1.set_ylim(0, clip_threshold)
        ax1.text(0.98, 0.98, f'Note: Y-axis clipped at {clip_threshold:.2f}s', 
                transform=ax1.transAxes, ha='right', va='top', 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Log scale
    for algo in ALGORITHMS:
        if plot_data[algo]['bits']:
            ax2.plot(plot_data[algo]['bits'], plot_data[algo]['times'], 
                    marker='o', label=f'{algo} (bonus)', color=colors[algo], linewidth=2, markersize=8)
    
    ax2.set_xlabel('Bit Length', fontsize=12)
    ax2.set_ylabel('Average Time (seconds, log scale)', fontsize=12)
    ax2.set_title(f'Bonus Implementation Performance - Log Scale', fontsize=14, fontweight='bold')
    ax2.set_yscale('log')
    ax2.legend(loc='upper left')
    ax2.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout() # type: ignore
    output_file = log_dir / f'bonus_performance_{min(bit_lengths)}_{max(bit_lengths)}bit.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight') # type: ignore
    plt.close() # type: ignore
    
    print(f"  📊 Bonus comparison: {output_file}")

if __name__ == "__main__":
    main()
