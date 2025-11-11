#!/usr/bin/env python3
"""
Quick graph generator - Run comparisons first, then use this to generate graphs
This avoids re-running long tests just to get graphs.
"""

import sys
import subprocess
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')
except ImportError:
    print("Error: matplotlib not installed. Run: pip3 install matplotlib")
    sys.exit(1)

def test_algorithm_quick(algo, bits):
    """Quick test - just first case to get average time."""
    test_file = Path(f'test_cases/{bits:02d}bit/case_1.txt')
    if not test_file.exists():
        return None
    
    script = Path(algo) / 'main_optimized.py'
    if not script.exists():
        return None
    
    try:
        import time
        start = time.time()
        result = subprocess.run(
            ['python3', str(script), str(test_file)],
            capture_output=True,
            text=True,
            timeout=10
        )
        elapsed = time.time() - start
        
        if result.returncode == 0 and ("Solution" in result.stdout or "PASSED" in result.stdout):
            return elapsed
        elif "requires exact point order" in result.stdout:
            return -1  # N/A
        return None
    except:
        return None

def generate_comparison_graphs(bit_start, bit_end):
    """Generate performance graphs for main algorithms."""
    
    ALGORITHMS = ['BruteForce', 'BabyStep', 'PohligHellman', 'PollardRho', 'LasVegas']
    colors = {
        'BruteForce': '#e74c3c',
        'BabyStep': '#3498db',
        'PohligHellman': '#2ecc71',
        'PollardRho': '#f39c12',
        'LasVegas': '#9b59b6'
    }
    
    print(f"Collecting performance data for {bit_start}-{bit_end} bits...")
    plot_data = {algo: {'bits': [], 'times': []} for algo in ALGORITHMS}
    
    for bits in range(bit_start, bit_end + 1):
        print(f"  Testing {bits}-bit...", end='', flush=True)
        for algo in ALGORITHMS:
            if bits > 18 and algo in ['PollardRho', 'LasVegas']:
                continue
            
            time_val = test_algorithm_quick(algo, bits)
            if time_val and time_val > 0:
                plot_data[algo]['bits'].append(bits)
                plot_data[algo]['times'].append(time_val)
        print(" ✓")
    
    # Create output directory
    output_dir = Path('graphs')
    output_dir.mkdir(exist_ok=True)
    
    # Generate plots
    print("\nGenerating graphs...")
    
    # Plot 1: Combined view
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Linear scale
    max_time = 0
    for algo in ALGORITHMS:
        if plot_data[algo]['bits']:
            times = plot_data[algo]['times']
            max_time = max(max_time, max(times) if times else 0)
            ax1.plot(plot_data[algo]['bits'], times, 
                    marker='o', label=algo, color=colors[algo], linewidth=2, markersize=6)
    
    ax1.set_xlabel('Bit Length', fontsize=12)
    ax1.set_ylabel('Time (seconds)', fontsize=12)
    ax1.set_title(f'ECC ECDLP Performance ({bit_start}-{bit_end} bits)', fontsize=14, fontweight='bold')
    ax1.legend(loc='upper left', fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Clip extreme values
    if max_time > 5:
        clip = max_time * 0.6
        ax1.set_ylim(0, clip)
        ax1.text(0.98, 0.98, f'Note: Clipped at {clip:.2f}s\n(removes outliers)', 
                transform=ax1.transAxes, ha='right', va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
    
    # Log scale
    for algo in ALGORITHMS:
        if plot_data[algo]['bits']:
            ax2.plot(plot_data[algo]['bits'], plot_data[algo]['times'], 
                    marker='o', label=algo, color=colors[algo], linewidth=2, markersize=6)
    
    ax2.set_xlabel('Bit Length', fontsize=12)
    ax2.set_ylabel('Time (seconds, log scale)', fontsize=12)
    ax2.set_title(f'ECC ECDLP Performance - Log Scale', fontsize=14, fontweight='bold')
    ax2.set_yscale('log')
    ax2.legend(loc='upper left', fontsize=10)
    ax2.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    output_file = output_dir / f'performance_{bit_start}_{bit_end}bit.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved: {output_file}")
    
    # Plot 2: Individual algorithms
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    
    for idx, algo in enumerate(ALGORITHMS):
        if plot_data[algo]['bits']:
            axes[idx].plot(plot_data[algo]['bits'], plot_data[algo]['times'], 
                          marker='o', color=colors[algo], linewidth=2.5, markersize=8)
            axes[idx].set_title(f'{algo}', fontsize=13, fontweight='bold', color=colors[algo])
            axes[idx].set_xlabel('Bit Length', fontsize=10)
            axes[idx].set_ylabel('Time (s)', fontsize=10)
            axes[idx].grid(True, alpha=0.3)
            
            if plot_data[algo]['times']:
                min_t = min(plot_data[algo]['times'])
                max_t = max(plot_data[algo]['times'])
                avg_t = sum(plot_data[algo]['times']) / len(plot_data[algo]['times'])
                axes[idx].text(0.05, 0.95, 
                              f'Min: {min_t:.4f}s\nAvg: {avg_t:.4f}s\nMax: {max_t:.4f}s', 
                              transform=axes[idx].transAxes, va='top', fontsize=9,
                              bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    
    axes[5].axis('off')
    
    plt.suptitle(f'Individual Algorithm Performance ({bit_start}-{bit_end} bits)', 
                fontsize=16, fontweight='bold')
    plt.tight_layout()
    output_file = output_dir / f'individual_{bit_start}_{bit_end}bit.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved: {output_file}")
    print(f"\n📊 Graphs saved in: {output_dir}/")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 generate_graphs.py <bit_start> [bit_end]")
        print("Example: python3 generate_graphs.py 10 18")
        sys.exit(1)
    
    bit_start = int(sys.argv[1])
    bit_end = int(sys.argv[2]) if len(sys.argv) > 2 else bit_start
    
    print("="*70)
    print("ECC ECDLP Performance Graph Generator")
    print("="*70)
    
    generate_comparison_graphs(bit_start, bit_end)
    
    print("\n" + "="*70)
    print("✓ Graph generation complete!")
    print("="*70)
