#!/usr/bin/env python3
"""Comparison script with visualization - tests algorithms and generates performance graphs"""

import subprocess
import time
import sys
from pathlib import Path

ALGORITHMS = ['BruteForce', 'BabyStep', 'PohligHellman'] # 'PollardRho', 'LasVegas']

# Try to import matplotlib, but don't fail if not available
try:
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Note: matplotlib not available, skipping graph generation")

def test_algorithm(algo, bit_length):
    """Test one algorithm on all 5 test cases for a bit length."""
    results = []
    
    for case_num in range(1, 6):
        test_file = Path(f'test_cases/{bit_length:02d}bit/case_{case_num}.txt')
        if not test_file.exists():
            results.append((None, None))
            continue
        
        script = Path(algo) / 'main_optimized.py'
        if not script.exists():
            results.append((None, None))
            continue
        
        try:
            start = time.time()
            result = subprocess.run(
                ['python3', str(script), str(test_file)],
                capture_output=True,
                text=True,
                timeout=15  # 15 seconds timeout
            )
            elapsed = time.time() - start
            
            # Extract attempt info for probabilistic algorithms
            attempts = None
            if algo in ['PollardRho', 'LasVegas']:
                if 'attempt' in result.stdout.lower() or 'Attempt' in result.stdout:
                    # Extract attempt number from output
                    import re
                    # Try different patterns
                    match = re.search(r'found in attempt (\d+)', result.stdout)
                    if not match:
                        match = re.search(r'Attempts?: (\d+)', result.stdout, re.IGNORECASE)
                    if match:
                        attempts = match.group(1)
            
            if result.returncode == 0 and ("Solution" in result.stdout or "PASSED" in result.stdout):
                results.append((elapsed, attempts))
            elif "requires exact point order" in result.stdout:
                results.append((-1, None))  # Special marker for N/A
            else:
                results.append((None, None))
        except subprocess.TimeoutExpired:
            results.append((None, None))
        except Exception as e:
            results.append((None, None))
    
    return results

def format_results(results, algo_name):
    """Format test results showing all 5 cases and average."""
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
        passed = f"({len(valid_times)}/5)"
    else:
        avg_str = "FAILED"
        passed = "(0/5)"
    
    # Add attempt info for probabilistic algorithms
    if attempts_list and algo_name in ['PollardRho', 'LasVegas']:
        avg_attempts = sum(int(a) for a in attempts_list) / len(attempts_list)
        avg_str += f" [~{avg_attempts:.0f} attempts]"
    
    return f"{time_str} | {avg_str:25s} {passed}"

def main():
    import sys
    
    # Parse arguments
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
    print("Format: Case1 | Case2 | Case3 | Case4 | Case5 | Average (passed/total)")
    print("Note: PohligHellman shows N/A when order factorization fails (expected for random curves)")
    print("      PollardRho/LasVegas show attempt count (probabilistic algorithms)")
    print("=" * 120)
    
    # Store data for plotting
    plot_data = {algo: {'bits': [], 'times': []} for algo in ALGORITHMS}
    
    for bits in range(bit_start, bit_end + 1):
        print(f"\n{bits}-bit:")
        for algo in ALGORITHMS:
            # Skip probabilistic algos for bits > 18
            if bits > 18 and algo in ['PollardRho', 'LasVegas']:
                print(f"  {algo:15s}: INFEASIBLE (not feasible beyond 18 bits)")
                continue
            
            results = test_algorithm(algo, bits)
            formatted = format_results(results, algo)
            print(f"  {algo:15s}: {formatted}")
            
            # Collect data for plotting
            valid_times = [t for t, _ in results if t is not None and t > 0]
            if valid_times:
                avg_time = sum(valid_times) / len(valid_times)
                plot_data[algo]['bits'].append(bits)
                plot_data[algo]['times'].append(avg_time)
    
    # Generate plots
    if HAS_MATPLOTLIB:
        print("\n" + "=" * 120)
        print("Generating performance graphs...")
        generate_plots(plot_data, bit_start, bit_end)
        print("✓ Graphs saved!")
    
    print("=" * 120)

def generate_plots(plot_data, bit_start, bit_end):
    """Generate performance visualization graphs."""
    
    # Create output directory
    output_dir = Path('graphs')
    output_dir.mkdir(exist_ok=True)
    
    # Color scheme for algorithms
    colors = {
        'BruteForce': '#e74c3c',      # Red
        'BabyStep': '#3498db',         # Blue
        'PohligHellman': '#2ecc71',    # Green
        'PollardRho': '#f39c12',       # Orange
        'LasVegas': '#9b59b6'          # Purple
    }
    
    # Plot 1: Linear scale (clipped for readability)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Linear scale plot
    max_time = 0
    for algo in ALGORITHMS:
        if plot_data[algo]['bits']:
            times = plot_data[algo]['times']
            max_time = max(max_time, max(times))
            ax1.plot(plot_data[algo]['bits'], times, 
                    marker='o', label=algo, color=colors[algo], linewidth=2)
    
    ax1.set_xlabel('Bit Length', fontsize=12)
    ax1.set_ylabel('Average Time (seconds)', fontsize=12)
    ax1.set_title(f'ECC ECDLP Performance ({bit_start}-{bit_end} bits) - Linear Scale', fontsize=14, fontweight='bold')
    ax1.legend(loc='upper left')
    ax1.grid(True, alpha=0.3)
    
    # Clip exceptionally high values for better visualization
    clip_threshold = max_time * 0.8 if max_time > 5 else None
    if clip_threshold:
        ax1.set_ylim(0, clip_threshold)
        ax1.text(0.98, 0.98, f'Note: Y-axis clipped at {clip_threshold:.2f}s', 
                transform=ax1.transAxes, ha='right', va='top', 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Log scale plot
    for algo in ALGORITHMS:
        if plot_data[algo]['bits']:
            ax2.plot(plot_data[algo]['bits'], plot_data[algo]['times'], 
                    marker='o', label=algo, color=colors[algo], linewidth=2)
    
    ax2.set_xlabel('Bit Length', fontsize=12)
    ax2.set_ylabel('Average Time (seconds, log scale)', fontsize=12)
    ax2.set_title(f'ECC ECDLP Performance ({bit_start}-{bit_end} bits) - Log Scale', fontsize=14, fontweight='bold')
    ax2.set_yscale('log')
    ax2.legend(loc='upper left')
    ax2.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    output_file = output_dir / f'comparison_{bit_start}_{bit_end}bit.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  📊 Main comparison: {output_file}")
    
    # Plot 2: Individual algorithm performance (separate subplots)
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    
    for idx, algo in enumerate(ALGORITHMS):
        if plot_data[algo]['bits']:
            axes[idx].plot(plot_data[algo]['bits'], plot_data[algo]['times'], 
                          marker='o', color=colors[algo], linewidth=2, markersize=8)
            axes[idx].set_title(f'{algo}', fontsize=12, fontweight='bold')
            axes[idx].set_xlabel('Bit Length', fontsize=10)
            axes[idx].set_ylabel('Avg Time (s)', fontsize=10)
            axes[idx].grid(True, alpha=0.3)
            
            # Add statistics
            if plot_data[algo]['times']:
                min_time = min(plot_data[algo]['times'])
                max_time = max(plot_data[algo]['times'])
                axes[idx].text(0.05, 0.95, 
                              f'Min: {min_time:.4f}s\nMax: {max_time:.4f}s', 
                              transform=axes[idx].transAxes, 
                              verticalalignment='top',
                              bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Hide the 6th subplot (we only have 5 algorithms)
    axes[5].axis('off')
    
    plt.suptitle(f'Individual Algorithm Performance ({bit_start}-{bit_end} bits)', 
                fontsize=16, fontweight='bold', y=1.00)
    plt.tight_layout()
    output_file = output_dir / f'individual_{bit_start}_{bit_end}bit.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  📊 Individual plots: {output_file}")

if __name__ == "__main__":
    main()
