#!/usr/bin/env python3
"""Comparison script with visualization - tests algorithms and generates performance graphs"""

import subprocess
import time
import sys
from pathlib import Path

# Enable the main algorithms for comparison
ALGORITHMS = ['BruteForce', 'BabyStep', 'PohligHellman', 'PollardRho', 'LasVegas']
#ALGORITHMS = ['PohligHellman','PollardRho', 'LasVegas']
#ALGORITHMS = ['BabyStep', 'PohligHellman', 'PollardRho', 'LasVegas']


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

# Try to import matplotlib
try:
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Note: matplotlib not available, skipping graph generation")

def test_algorithm(algo, bit_length):
    """Test one algorithm on all available test cases for a bit length."""
    results = []
    case_files = discover_case_files(bit_length)
    
    # Process files
    for test_file in case_files:
        script = Path(algo) / 'main_optimized.py'
        if not script.exists():
            # Try main.py if optimized doesn't exist
            script = Path(algo) / 'main.py'
            if not script.exists():
                results.append((None, None))
                continue
        
        try:
            start = time.time()
            result = subprocess.run(
                [sys.executable, str(script), str(test_file)], # Use sys.executable for venv safety
                capture_output=True,
                text=True,
                timeout=60  # 60 seconds timeout per case
            )
            elapsed = time.time() - start
            
            # Extract attempt info for PollardRho
            attempts = None
            if algo in ['PollardRho', 'LasVegas', 'MOV']:
                import re
                match = re.search(r'attempt (\d+)', result.stdout, re.IGNORECASE)
                if match:
                    attempts = match.group(1)
            
            if result.returncode == 0 and ("Solution" in result.stdout or "PASSED" in result.stdout):
                results.append((elapsed, attempts))
            else:
                results.append((None, None)) # Failed or wrong answer
                
        except subprocess.TimeoutExpired:
            results.append((None, None)) # Timeout
        except Exception:
            results.append((None, None)) # Crash
    
    return results

def format_results(results, algo_name):
    """Format test results showing all cases and average."""
    valid_times = [r[0] for r in results if r[0] is not None and r[0] > 0]
    attempts_list = [r[1] for r in results if r[1] is not None]
    
    # Format individual times
    time_parts = []
    for r in results:
        time_val, attempt_val = r
        if time_val and time_val > 0:
            time_str = f"{time_val*1000:5.1f}ms"
        elif time_val is None:
            time_str = "FAIL "
        else:
            time_str = " N/A "
        time_parts.append(time_str)
    
    time_str = " | ".join(time_parts)
    
    # Calculate average
    if valid_times:
        avg = sum(valid_times) / len(valid_times)
        if avg < 1:
            avg_str = f"avg={avg*1000:.1f}ms"
        else:
            avg_str = f"avg={avg:.3f}s"
        denom = len(results) if results else 0
        passed = f"({len(valid_times)}/{denom})"
    else:
        avg_str = "FAILED"
        denom = len(results) if results else 0
        passed = f"(0/{denom})"
    
    # Add attempt info
    if attempts_list and algo_name in ['PollardRho', 'LasVegas']:
        try:
            avg_attempts = sum(int(a) for a in attempts_list) / len(attempts_list)
            avg_str += f" [~{avg_attempts:.1f} tries]"
        except:
            pass
    
    return f"{time_str} | {avg_str:25s} {passed}"

def main():
    import sys
    
    if len(sys.argv) > 1:
        try:
            bit_start = int(sys.argv[1])
            bit_end = int(sys.argv[2]) if len(sys.argv) > 2 else bit_start
        except:
            print("Usage: python3 run_comparison.py [bit_start] [bit_end]")
            sys.exit(1)
    else:
        bit_start, bit_end = 10, 30
    
    print("=" * 120)
    print(f"ECC ECDLP Performance Comparison ({bit_start}-{bit_end} bits)")
    print("=" * 120)
    
    plot_data = {algo: {'bits': [], 'times': []} for algo in ALGORITHMS}
    
    for bits in range(bit_start, bit_end + 1):
        print(f"\n{bits}-bit:")
        for algo in ALGORITHMS:
            
            # --- SMART LIMITS ---
            if algo == 'BruteForce' and bits > 24:
                print(f"  {algo:15s}: SKIPPED (exponential time)")
                continue
            if algo == 'BabyStep' and bits > 50:
                print(f"  {algo:15s}: SKIPPED (memory limit)")
                continue
            # PollardRho and PohligHellman allowed on all sizes
            
            results = test_algorithm(algo, bits)
            formatted = format_results(results, algo)
            print(f"  {algo:15s}: {formatted}")
            
            # Collect data
            valid_times = [t for t, _ in results if t is not None and t > 0]
            if valid_times:
                avg_time = sum(valid_times) / len(valid_times)
                plot_data[algo]['bits'].append(bits)
                plot_data[algo]['times'].append(avg_time)
    
    # Generate plots
    if HAS_MATPLOTLIB:
        print("\n" + "=" * 120)
        generate_plots(plot_data, bit_start, bit_end)
        print("✓ Graphs saved to 'graphs/' folder")
    
    print("=" * 120)

def generate_plots(plot_data, bit_start, bit_end):
    """Generate comprehensive performance graphs."""
    output_dir = Path('graphs')
    output_dir.mkdir(exist_ok=True)
    
    # Standard Colors
    colors = {
        'BruteForce': '#e74c3c',      # Red
        'BabyStep': '#3498db',         # Blue
        'PohligHellman': '#2ecc71',    # Green
        'PollardRho': '#f39c12',       # Orange
        'LasVegas': '#9b59b6',         # Purple
        'MOV': '#34495e'               # Dark Grey
    }
    
    # ---------------------------------------------------------
    # GRAPH 1: Comparison (Linear & Log Side-by-Side)
    # ---------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Linear Plot
    max_time = 0
    for algo in ALGORITHMS:
        if plot_data[algo]['bits']:
            times = plot_data[algo]['times']
            if times:
                max_time = max(max_time, max(times))
            ax1.plot(plot_data[algo]['bits'], times, 
                    marker='o', label=algo, color=colors.get(algo, 'black'), linewidth=2)
    
    ax1.set_xlabel('Bit Length')
    ax1.set_ylabel('Time (s)')
    ax1.set_title('Linear Scale')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Clip linear graph if one algo explodes (e.g. BruteForce)
    if max_time > 10: 
        ax1.set_ylim(0, 10)
        ax1.text(0.02, 0.98, 'Clipped at 10s', transform=ax1.transAxes, 
                 bbox=dict(facecolor='yellow', alpha=0.5), va='top')

    # Log Plot
    for algo in ALGORITHMS:
        if plot_data[algo]['bits']:
            ax2.plot(plot_data[algo]['bits'], plot_data[algo]['times'], 
                    marker='o', label=algo, color=colors.get(algo, 'black'), linewidth=2)
            
    ax2.set_xlabel('Bit Length')
    ax2.set_ylabel('Time (s)')
    ax2.set_title('Log Scale')
    ax2.set_yscale('log')
    ax2.legend()
    ax2.grid(True, which="both", alpha=0.2)
    
    plt.tight_layout()
    plt.savefig(output_dir / f'comparison_{bit_start}_{bit_end}bit.png')
    plt.close()
    
    # ---------------------------------------------------------
    # GRAPH 2: Individual Algorithms
    # ---------------------------------------------------------
    num_algos = len(ALGORITHMS)
    cols = 3
    rows = (num_algos + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
    axes = axes.flatten()
    
    for idx, algo in enumerate(ALGORITHMS):
        if idx < len(axes):
            ax = axes[idx]
            if plot_data[algo]['bits']:
                ax.plot(plot_data[algo]['bits'], plot_data[algo]['times'], 
                       marker='o', color=colors.get(algo, 'black'), linewidth=2)
                ax.set_title(algo, fontweight='bold')
                ax.set_xlabel('Bit Length')
                ax.set_ylabel('Time (s)')
                ax.grid(True, alpha=0.3)
            else:
                ax.text(0.5, 0.5, "No Data", ha='center', va='center')
                ax.set_title(algo)
                
    # Hide unused subplots
    for i in range(len(ALGORITHMS), len(axes)):
        axes[i].axis('off')
        
    plt.tight_layout()
    plt.savefig(output_dir / f'individual_{bit_start}_{bit_end}bit.png')
    plt.close()

if __name__ == "__main__":
    main()