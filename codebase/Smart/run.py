#!/usr/bin/env python3
import subprocess
import time
import sys
import argparse
from pathlib import Path
import matplotlib.pyplot as plt 

# Map curve types to the algorithm that breaks them
SCENARIOS = [

    {
        "name": "Anomalous (Smart)",
        "folder": "Anomalous",
        "script": "smart.py",
        "color": "orange",
        "marker": "^"
    },
    {
        "name": "Generic (Pollard Rho)",
        "folder": "Generic",
        "script": "../PollardRho/main_optimized.py",
        "color": "blue",
        "marker": "s"
    },
    {
        "name": "Smooth (Pohlig-Hellman)",
        "folder": "PH_friendly",
        "script": "../PohligHellman/main_optimized.py",
        "color": "green",
        "marker": "x"
    }
]

def run_script(script_path, input_file):
    start = time.time()
    try:
        # 10 second timeout for attacks
        subprocess.run(
            ['python3', script_path, str(input_file)],
            capture_output=True, text=True, timeout=10
        )
        return time.time() - start
    except subprocess.TimeoutExpired:
        return 10.0  # Cap at timeout
    except Exception:
        return 10.0

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="ECC Security Benchmark")
    parser.add_argument("--start", type=int, default=10, help="Starting bit size (default: 10)")
    parser.add_argument("--end", type=int, default=20, help="Ending bit size (default: 20)")
    parser.add_argument("--step", type=int, default=1, help="Step size for bits (default: 1)")
    args = parser.parse_args()
    
    START_BIT = args.start
    END_BIT = args.end
    STEP = args.step
    
    base_path = Path(__file__).parent
    test_cases_dir = base_path / "test_cases"
    
    # Check if test cases exist
    if not test_cases_dir.exists():
        print("[ERROR] Test cases not found!")
        print("Please run: python3 Smart/gen.py --start 20 --end 40 --step 2")
        sys.exit(1)
    
    print("=" * 80)
    print("ECC SECURITY BENCHMARK: Vulnerable vs. Secure Curves")
    print("=" * 80)
    print(f"Bit range: {START_BIT} to {END_BIT} | Testing 3 attack types on 3 curve types")
    print("=" * 80)
    print()
    
    # 1. Run Benchmarks
    # Structure: results[curve_type][attack_name] = {'bits': [], 'times': []}
    curve_types = ["Anomalous", "Generic", "PH_friendly"]
    attacks = {
        "Smart": "smart.py",
        "Pollard Rho": "../PollardRho/main_optimized.py",
        "Pohlig-Hellman": "../PohligHellman/main_optimized.py"
    }
    
    # Initialize results structure
    results_cross = {}
    for curve in curve_types:
        results_cross[curve] = {}
        for attack in attacks:
            results_cross[curve][attack] = {'bits': [], 'times': []}
    
    # Also keep old structure for display
    results = {s['name']: {'bits': [], 'times': []} for s in SCENARIOS}
    
    print("Running benchmarks...")
    print()
    
    for bits in range(START_BIT, END_BIT + 1, STEP):
        print(f"[{bits}-bit curves]")
        
        # First pass: optimal attacks (for display)
        for sc in SCENARIOS:
            case_dir = test_cases_dir / sc['folder'] / f"{bits}bit"
            
            if not case_dir.exists():
                print(f"  {sc['name']:<30} SKIP (no test cases)")
                continue
            
            cases = sorted(list(case_dir.glob("case_*.txt")))
            
            if not cases:
                print(f"  {sc['name']:<30} SKIP (no test cases)")
                continue
            
            times = []
            script_full_path = str(base_path / sc['script'])
            
            for case_file in cases:
                time_taken = run_script(script_full_path, case_file)
                times.append(time_taken)
            
            avg_time = sum(times) / len(times) if times else 10.0
            min_time = min(times) if times else 10.0
            max_time = max(times) if times else 10.0
            
            results[sc['name']]['bits'].append(bits)
            results[sc['name']]['times'].append(avg_time)
            
            print(f"  {sc['name']:<30} {avg_time:7.4f}s (n={len(times)}, min={min_time:.4f}s, max={max_time:.4f}s)")
        
        # Second pass: cross-testing ALL attacks on ALL curve types
        for curve_type in curve_types:
            case_dir = test_cases_dir / curve_type / f"{bits}bit"
            if not case_dir.exists():
                continue
            
            cases = sorted(list(case_dir.glob("case_*.txt")))
            if not cases:
                continue
            
            for attack_name, attack_script in attacks.items():
                times = []
                script_full_path = str(base_path / attack_script)
                
                for case_file in cases:
                    time_taken = run_script(script_full_path, case_file)
                    times.append(time_taken)
                
                avg_time = sum(times) / len(times) if times else 10.0
                results_cross[curve_type][attack_name]['bits'].append(bits)
                results_cross[curve_type][attack_name]['times'].append(avg_time)
        
        print()

    # 2. Plotting
    print()
    print("Generating graphs...")
    
    # Graph 1: Overall comparison (all attacks on one plot)
    plt.figure(figsize=(12, 7))
    
    for sc in SCENARIOS:
        data = results[sc['name']]
        if data['bits']:
            plt.plot(
                data['bits'], 
                data['times'], 
                label=sc['name'],
                color=sc['color'],
                marker=sc['marker'],
                linewidth=2,
                markersize=6
            )

    plt.xlabel("Key Size (Bits)", fontsize=12)
    plt.ylabel("Time to Break (Seconds)", fontsize=12)
    plt.title("ECC Security: All Attack Types Comparison", fontsize=14, fontweight='bold')
    plt.axhline(y=10.0, color='gray', linestyle='--', alpha=0.5, label="Timeout (10s)")
    
    plt.yscale('log')
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend(fontsize=10)
    plt.tight_layout()
    
    out_file_all = base_path / "benchmark_all_attacks.png"
    plt.savefig(out_file_all, dpi=300)
    plt.close()
    
    curve_labels = {
        "Anomalous": "Anomalous",
        "Generic": "Generic (Secure)",
        "PH_friendly": "Smooth Order"
    }
    curve_colors = {
        "Anomalous": "orange", 
        "Generic": "blue",
        "PH_friendly": "green"
    }
    
    # Graph 2: Pollard Rho on ALL Curve Types
    plt.figure(figsize=(12, 7))
    
    for curve_type in curve_types:
        data = results_cross[curve_type].get("Pollard Rho", {'bits': [], 'times': []})
        if data['bits']:
            lw = 3 if curve_type == "Generic" else 2
            alpha = 1.0 if curve_type == "Generic" else 0.6
            ls = '-' if curve_type == "Generic" else '--'
            plt.plot(data['bits'], data['times'], 
                    label=curve_labels[curve_type].replace("(Secure)", "(Works on all)"),
                    color=curve_colors[curve_type],
                    marker='s', linewidth=lw, markersize=6, alpha=alpha, linestyle=ls)
    
    plt.xlabel("Key Size (Bits)", fontsize=12)
    plt.ylabel("Time to Break (Seconds)", fontsize=12)
    plt.title("Pollard Rho Performance: Universal Attack on All Curve Types", fontsize=14, fontweight='bold')
    plt.yscale('log')
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend(fontsize=10)
    plt.tight_layout()
    
    out_file_rho = base_path / "benchmark_pollard_comparison.png"
    plt.savefig(out_file_rho, dpi=300)
    plt.close()
    
    # Graph 3: Pohlig-Hellman on ALL Curve Types
    plt.figure(figsize=(12, 7))
    
    for curve_type in curve_types:
        data = results_cross[curve_type].get("Pohlig-Hellman", {'bits': [], 'times': []})
        if data['bits']:
            lw = 3 if curve_type == "PH_friendly" else 2
            alpha = 1.0 if curve_type == "PH_friendly" else 0.6
            ls = '-' if curve_type == "PH_friendly" else '--'
            plt.plot(data['bits'], data['times'], 
                    label=curve_labels[curve_type].replace("Smooth Order", "Smooth Order (Vulnerable to PH)"),
                    color=curve_colors[curve_type],
                    marker='x', linewidth=lw, markersize=8, alpha=alpha, linestyle=ls)
    
    plt.xlabel("Key Size (Bits)", fontsize=12)
    plt.ylabel("Time to Break (Seconds)", fontsize=12)
    plt.title("Pohlig-Hellman Performance: Testing on All Curve Types", fontsize=14, fontweight='bold')
    plt.yscale('log')
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend(fontsize=10)
    plt.tight_layout()
    
    out_file_ph = base_path / "benchmark_pohlig_comparison.png"
    plt.savefig(out_file_ph, dpi=300)
    plt.close()
    
    print()
    print("=" * 80)
    print("✓ Benchmark complete! Generated 3 comparison graphs:")
    print(f"  • Overall:        {out_file_all}")
    print(f"  • Pollard Focus:  {out_file_rho}")
    print(f"  • Pohlig Focus:   {out_file_ph}")
    print("=" * 80)

if __name__ == "__main__":
    main()