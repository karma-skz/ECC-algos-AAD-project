#!/usr/bin/env python3
import subprocess
import time
import sys
from pathlib import Path
import matplotlib.pyplot as plt

# --- CONFIGURATION ---
START_BIT = 20
END_BIT = 40
STEP = 2  # Step size for bits (20, 22, 24...)

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
    base_path = Path(__file__).parent
    test_cases_dir = base_path / "test_cases"
    
    # 1. Generate Data
    print(f"[-] Generating curves from {START_BIT} to {END_BIT} bits...")
    subprocess.run([
        "python3", "gen.py", 
        "--start", str(START_BIT), 
        "--end", str(END_BIT), 
        "--step", str(STEP)
    ], cwd=base_path)

    # 2. Run Benchmarks
    results = {s['name']: {'bits': [], 'times': []} for s in SCENARIOS}
    
    print("\n[-] Running benchmarks...")
    for bits in range(START_BIT, END_BIT + 1, STEP):
        print(f"  Testing {bits}-bit curves...")
        
        for sc in SCENARIOS:
            # Find a test case for this specific bit & type
            case_dir = test_cases_dir / sc['folder'] / f"{bits}bit"
            cases = list(case_dir.glob("case_*.txt"))
            
            if not cases:
                continue
                
            # Run test on the first available case
            script_full_path = str(base_path / sc['script'])
            time_taken = run_script(script_full_path, cases[0])
            
            results[sc['name']]['bits'].append(bits)
            results[sc['name']]['times'].append(time_taken)
            print(f"    {sc['name']:<25}: {time_taken:.4f}s")

    # 3. Plotting
    print("\n[-] Generating Graph...")
    plt.figure(figsize=(10, 6))
    
    for sc in SCENARIOS:
        data = results[sc['name']]
        if data['bits']:
            plt.plot(
                data['bits'], 
                data['times'], 
                label=sc['name'],
                color=sc['color'],
                marker=sc['marker'],
                linewidth=2
            )

    plt.xlabel("Key Size (Bits)", fontsize=12)
    plt.ylabel("Time to Break (Seconds)", fontsize=12)
    plt.title("ECC Security Gap: Weak vs. Strong Curves", fontsize=14, fontweight='bold')
    plt.axhline(y=10.0, color='gray', linestyle='--', alpha=0.5, label="Timeout (10s)")
    
    plt.yscale('log')  # Important for visualizing exponential diff
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend()
    
    out_file = base_path / "benchmark_results.png"
    plt.savefig(out_file, dpi=300)
    print(f"[Image of line graph] Graph saved to: {out_file}")

if __name__ == "__main__":
    main()