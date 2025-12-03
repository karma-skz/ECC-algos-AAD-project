#!/usr/bin/env python3
"""run_comparison.py
Compares ECC ECDLP algorithms across test cases and generates graphs.

Usage:
    python3 run_comparison.py [bit_start] [bit_end]

Produces:
    - graphs/comparison_<start>_<end>bit.png
    - graphs/individual_<start>_<end>bit.png
    - graphs/<ALGO>_cases_<start>_<end>bit.png  (per-case Linear+Log for each algorithm)
"""
import subprocess
import time
import sys
from pathlib import Path
import re

# Algorithms to test
ALGORITHMS = ['BruteForce', 'BabyStep', 'PohligHellman', 'PollardRho', 'LasVegas']

# Try to import matplotlib
try:
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')
    HAS_MATPLOTLIB = True
except Exception:
    HAS_MATPLOTLIB = False
    print("Note: matplotlib not available, skipping graph generation")

# -------------------------
# Utilities
# -------------------------
def discover_case_files(bit_length: int):
    """Return sorted list of case files for given bitsize."""
    cases_dir = Path('test_cases') / f'{bit_length:02d}bit'
    if not cases_dir.exists():
        return []
    def case_key(p: Path):
        m = re.search(r'case_(\d+)\.txt$', p.name)
        return int(m.group(1)) if m else 0
    return sorted(cases_dir.glob('case_*.txt'), key=case_key)

def test_algorithm(algo, bit_length):
    """Test one algorithm on all available test cases for a bit length.
    Returns list of (elapsed_seconds or None, attempts or None) for each case file.
    """
    results = []
    case_files = discover_case_files(bit_length)
    for test_file in case_files:
        script = Path(algo) / 'main_optimized.py'
        if not script.exists():
            script = Path(algo) / 'main.py'
            if not script.exists():
                results.append((None, None))
                continue

        try:
            start = time.time()
            proc = subprocess.run(
                [sys.executable, str(script), str(test_file)],
                capture_output=True, text=True, timeout=60
            )
            elapsed = time.time() - start

            attempts = None
            if algo in ['PollardRho', 'LasVegas', 'MOV']:
                m = re.search(r'attempt(?:s)?\s*(\d+)', proc.stdout, re.IGNORECASE)
                if m:
                    attempts = m.group(1)

            # Heuristic success detection
            if proc.returncode == 0 and ("Solution" in proc.stdout or "PASSED" in proc.stdout or "d =" in proc.stdout):
                results.append((elapsed, attempts))
            else:
                results.append((None, None))
        except subprocess.TimeoutExpired:
            results.append((None, None))
        except Exception:
            results.append((None, None))
    return results

def format_results(results, algo_name):
    """Format test results showing all cases and average."""
    valid_times = [r[0] for r in results if r[0] is not None and r[0] > 0]
    attempts_list = [r[1] for r in results if r[1] is not None]

    time_parts = []
    for r in results:
        time_val, _ = r
        if time_val and time_val > 0:
            time_parts.append(f"{time_val*1000:5.1f}ms")
        elif time_val is None:
            time_parts.append("FAIL ")
        else:
            time_parts.append(" N/A ")

    time_str = " | ".join(time_parts)

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

    if attempts_list and algo_name in ['PollardRho', 'LasVegas']:
        try:
            avg_attempts = sum(int(a) for a in attempts_list) / len(attempts_list)
            avg_str += f" [~{avg_attempts:.1f} tries]"
        except Exception:
            pass

    return f"{time_str} | {avg_str:25s} {passed}"

# -------------------------
# Plot generation
# -------------------------
def generate_plots(plot_data, bit_start, bit_end):
    """Generate comparison and individual algorithm grid plots."""
    output_dir = Path('graphs')
    output_dir.mkdir(exist_ok=True)

    colors = {
        'BruteForce': '#e74c3c',
        'BabyStep': '#3498db',
        'PohligHellman': '#2ecc71',
        'PollardRho': '#f39c12',
        'LasVegas': '#9b59b6',
        'MOV': '#34495e'
    }

    # GRAPH 1: Comparison Linear + Log
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    max_time = 0
    for algo in ALGORITHMS:
        bits = plot_data[algo]['bits']
        times = plot_data[algo]['times']
        if bits and times:
            ax1.plot(bits, times, marker='o', label=algo, color=colors.get(algo, 'black'), linewidth=2)
            ax2.plot(bits, times, marker='o', label=algo, color=colors.get(algo, 'black'), linewidth=2)
            max_time = max(max_time, max(times))
    ax1.set_xlabel('Bit Length'); ax1.set_ylabel('Time (s)'); ax1.set_title('Linear Scale'); ax1.legend(); ax1.grid(True, alpha=0.3)
    if max_time > 10:
        ax1.set_ylim(0, 10)
        ax1.text(0.02, 0.98, 'Clipped at 10s', transform=ax1.transAxes, bbox=dict(facecolor='yellow', alpha=0.5), va='top')
    ax2.set_xlabel('Bit Length'); ax2.set_ylabel('Time (s)'); ax2.set_title('Log Scale'); ax2.set_yscale('log'); ax2.legend(); ax2.grid(True, which="both", alpha=0.2)
    plt.tight_layout()
    plt.savefig(output_dir / f'comparison_{bit_start}_{bit_end}bit.png')
    plt.close()

    # GRAPH 2: Individual algorithms grid
    num_algos = len(ALGORITHMS)
    cols = 3
    rows = (num_algos + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
    axes = axes.flatten()
    for idx, algo in enumerate(ALGORITHMS):
        ax = axes[idx]
        bits = plot_data[algo]['bits']
        times = plot_data[algo]['times']
        if bits and times:
            ax.plot(bits, times, marker='o', color=colors.get(algo, 'black'), linewidth=2)
            ax.set_title(algo, fontweight='bold')
            ax.set_xlabel('Bit Length'); ax.set_ylabel('Time (s)'); ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, "No Data", ha='center', va='center')
            ax.set_title(algo)
    for i in range(len(ALGORITHMS), len(axes)):
        axes[i].axis('off')
    plt.tight_layout()
    plt.savefig(output_dir / f'individual_{bit_start}_{bit_end}bit.png')
    plt.close()

def generate_case_plots(plot_data, bit_start, bit_end):
    """For each algorithm, produce a figure with Linear and Log subplots showing per-case curves."""
    output_dir = Path('graphs')
    output_dir.mkdir(exist_ok=True)

    for algo in ALGORITHMS:
        cases_dict = plot_data[algo].get('cases', {})
        bits = sorted(cases_dict.keys())
        if not bits:
            continue

        # Determine maximum number of cases across bits
        max_cases = max(len(cases_dict[b]) for b in bits)

        # Build series per case index
        case_series = []
        for case_idx in range(max_cases):
            series = []
            for b in bits:
                lst = cases_dict[b]
                if case_idx < len(lst):
                    series.append(lst[case_idx])
                else:
                    series.append(None)
            case_series.append(series)

        # Plot
        fig, (ax_lin, ax_log) = plt.subplots(1, 2, figsize=(16, 6))
        for idx, series in enumerate(case_series):
            ax_lin.plot(bits, series, marker='o', label=f'Case {idx+1}')
        ax_lin.set_title(f"{algo} — Linear Scale"); ax_lin.set_xlabel("Bit Length"); ax_lin.set_ylabel("Time (s)"); ax_lin.grid(True, alpha=0.3); ax_lin.legend()

        for idx, series in enumerate(case_series):
            ax_log.plot(bits, series, marker='o', label=f'Case {idx+1}')
        ax_log.set_title(f"{algo} — Log Scale"); ax_log.set_xlabel("Bit Length"); ax_log.set_yscale("log"); ax_log.grid(True, which="both", alpha=0.3); ax_log.legend()

        plt.tight_layout()
        plt.savefig(output_dir / f'{algo}_cases_{bit_start}_{bit_end}bit.png')
        plt.close()


def generate_casewise_all_algos_plots(plot_data, bit_start, bit_end):
    """For each case index, plot all algorithms on one graph across bit lengths."""
    output_dir = Path('graphs')
    output_dir.mkdir(exist_ok=True)

    # Determine total number of cases from ANY algorithm
    # (Take the maximum across all algos and bits)
    max_cases = 0
    for algo in ALGORITHMS:
        for b in plot_data[algo]['cases']:
            max_cases = max(max_cases, len(plot_data[algo]['cases'][b]))

    colors = {
        'BruteForce': '#e74c3c',
        'BabyStep': '#3498db',
        'PohligHellman': '#2ecc71',
        'PollardRho': '#f39c12',
        'LasVegas': '#9b59b6',
    }

    # For each case index generate a single plot with all algorithms
    for case_idx in range(max_cases):
        fig, (ax_lin, ax_log) = plt.subplots(1, 2, figsize=(16, 6))

        for algo in ALGORITHMS:
            bits = sorted(plot_data[algo]['cases'].keys())
            series = []

            for b in bits:
                lst = plot_data[algo]['cases'][b]
                if case_idx < len(lst):
                    series.append(lst[case_idx])
                else:
                    series.append(None)

            if any(v is not None for v in series):
                ax_lin.plot(bits, series, marker='o', label=algo, color=colors.get(algo, 'black'))
                ax_log.plot(bits, series, marker='o', label=algo, color=colors.get(algo, 'black'))

        # Linear plot settings
        ax_lin.set_title(f"Case {case_idx+1} — All Algorithms (Linear)")
        ax_lin.set_xlabel("Bit Length")
        ax_lin.set_ylabel("Time (s)")
        ax_lin.grid(True, alpha=0.3)
        ax_lin.legend()

        # Log plot settings
        ax_log.set_title(f"Case {case_idx+1} — All Algorithms (Log)")
        ax_log.set_xlabel("Bit Length")
        ax_log.set_yscale("log")
        ax_log.grid(True, which="both", alpha=0.3)
        ax_log.legend()

        plt.tight_layout()
        plt.savefig(output_dir / f'case_{case_idx+1}_algo_{bit_start}_{bit_end}bit.png')
        plt.close()

# -------------------------
# Main driver
# -------------------------
def main():
    if len(sys.argv) > 1:
        try:
            bit_start = int(sys.argv[1])
            bit_end = int(sys.argv[2]) if len(sys.argv) > 2 else bit_start
        except Exception:
            print("Usage: python3 run_comparison.py [bit_start] [bit_end]")
            sys.exit(1)
    else:
        bit_start, bit_end = 10, 30

    print("=" * 120)
    print(f"ECC ECDLP Performance Comparison ({bit_start}-{bit_end} bits)")
    print("=" * 120)

    # Initialize plot_data with per-algo buckets including per-case storage
    plot_data = {algo: {'bits': [], 'times': [], 'cases': {}} for algo in ALGORITHMS}

    for bits in range(bit_start, bit_end + 1):
        print(f"\n{bits}-bit:")
        for algo in ALGORITHMS:
            # SMART LIMITS
            if algo == 'BruteForce' and bits > 20:
                print(f"  {algo:15s}: SKIPPED (exponential time)")
                continue
            if algo == 'BabyStep' and bits > 50:
                print(f"  {algo:15s}: SKIPPED (memory limit)")
                continue

            results = test_algorithm(algo, bits)
            formatted = format_results(results, algo)
            print(f"  {algo:15s}: {formatted}")

            # Collect per-bit averages for summary plots
            valid_times = [t for t, _ in results if t is not None and t > 0]
            if valid_times:
                avg_time = sum(valid_times) / len(valid_times)
                plot_data[algo]['bits'].append(bits)
                plot_data[algo]['times'].append(avg_time)

            # Collect per-case times (preserve None for fails)
            case_times = [t for t, _ in results]
            # store only if there are any case files; else skip
            if case_times:
                plot_data[algo]['cases'][bits] = case_times

    # Generate plots if matplotlib present
    if HAS_MATPLOTLIB:
        print("\n" + "=" * 120)
        generate_plots(plot_data, bit_start, bit_end)
        generate_case_plots(plot_data, bit_start, bit_end)
        generate_casewise_all_algos_plots(plot_data, bit_start, bit_end)
        print("✓ Graphs saved to 'graphs/' folder")

    print("=" * 120)

if __name__ == "__main__":
    main()
