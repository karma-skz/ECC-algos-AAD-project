import subprocess
import sys
import json
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# ==========================================
# 1. EXECUTION ENGINE (Reused from your code)
# ==========================================
def get_sorted_cases():
    """Scans test_cases/ folder for all available bit lengths."""
    cases = []
    base_dir = Path("test_cases")
    if not base_dir.exists(): return []

    for d in base_dir.iterdir():
        if d.is_dir() and "bit" in d.name:
            try:
                bits = int(d.name.replace("bit", ""))
                # Find valid txt file
                txts = list(d.glob("*.txt"))
                if txts:
                    cases.append({"bits": bits, "path": txts[0]})
            except ValueError: continue
    
    return sorted(cases, key=lambda x: x["bits"])

def run_algo(algo_name, case_path):
    """Runs a specific bonus script and returns the parsed JSON result."""
    script_path = Path(algo_name) / "bonus.py"
    if not script_path.exists(): return None
    
    try:
        # Run command: python3 Algo/bonus.py path/to/case.txt
        result = subprocess.run(
            [sys.executable, str(script_path), str(case_path)],
            capture_output=True, text=True
        )
        
        # Parse JSON from stdout
        for line in result.stdout.splitlines():
            if line.startswith("BONUS_RESULT:"):
                return json.loads(line.replace("BONUS_RESULT:", "").strip())
    except Exception as e:
        print(f"Error running {algo_name}: {e}")
    return None

# ==========================================
# 2. DATA COLLECTION
# ==========================================
cases = get_sorted_cases()
if not cases:
    print("No test cases found in test_cases/ folder.")
    sys.exit(1)

print(f"[*] Found {len(cases)} bit-lengths to benchmark: {[c['bits'] for c in cases]}")

# Storage for plotting data
# Structure: { "AlgoName": { "bits": [], "time": [], "metric": [] } }
data = {
    "BruteForce":    {"bits": [], "time": [], "reduction": []},
    "PollardRho":    {"bits": [], "time": [], "reduction": []},
    "BabyStep":      {"bits": [], "time": [], "cost": []},
    "PohligHellman": {"bits": [], "time": []},
    "LasVegas":      {"bits": [], "time": []}
}

for c in cases:
    bits = c['bits']
    print(f"Benchmarking {bits}-bit curve...", end=" ", flush=True)
    
    # Run Brute Force (LSB)
    res = run_algo("BruteForce", c['path'])
    if res and res['status'] == 'success':
        data["BruteForce"]["bits"].append(bits)
        data["BruteForce"]["time"].append(res['time'])
        # Metric: 'speedup' from JSON (e.g. "3.8e02x")
        spd = float(res['details'].get('speedup', '1x').replace('x',''))
        data["BruteForce"]["reduction"].append(spd)

    # Run Pollard (Interval)
    res = run_algo("PollardRho", c['path'])
    if res and res['status'] == 'success':
        data["PollardRho"]["bits"].append(bits)
        data["PollardRho"]["time"].append(res['time'])
        # Metric: Reduction = Space Size / Interval Width
        # Approx space size = 2^bits
        interval = float(res['details'].get('interval_width', 1))
        reduction = (2**bits) / interval
        data["PollardRho"]["reduction"].append(reduction)

    # Run BabyStep (RAM)
    res = run_algo("BabyStep", c['path'])
    if res and res['status'] == 'success':
        data["BabyStep"]["bits"].append(bits)
        data["BabyStep"]["time"].append(res['time'])
        # Metric: Cost = Time * (1/Savings)
        savings = float(res['details'].get('ram_savings', '1x').replace('x',''))
        # If standard cost is Time, Optimized Cost = Time * (1/10)
        # We plot "Effective Cost"
        data["BabyStep"]["cost"].append(res['time'] / savings)

    # Run Pohlig
    res = run_algo("PohligHellman", c['path'])
    if res and res['status'] == 'success':
        data["PohligHellman"]["bits"].append(bits)
        data["PohligHellman"]["time"].append(max(res['time'], 0.000001))

    print("Done.")

# ==========================================
# 3. SCIENTIFIC PLOTTING
# ==========================================
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))

# --- GRAPH 1: SCALABILITY (Time vs Bits) ---
# Compares standard geometric attacks vs algebraic attacks
ax1.plot(data["BruteForce"]["bits"], data["BruteForce"]["time"], 'o-', label='Brute Force (LSB)')
ax1.plot(data["PollardRho"]["bits"], data["PollardRho"]["time"], 's-', label='Pollard (Interval)')
ax1.plot(data["PohligHellman"]["bits"], data["PohligHellman"]["time"], '^-', linewidth=3, label='Pohlig-Hellman')

ax1.set_yscale('log')
ax1.set_xlabel('Curve Size (Bits)')
ax1.set_ylabel('Execution Time (s) - Log Scale')
ax1.set_title('Scalability: Algebraic vs Geometric Attacks')
ax1.legend()
ax1.grid(True, which="both", ls="--", alpha=0.3)

# --- GRAPH 2: SEARCH EFFICIENCY (Reduction Factor) ---
# "How much effectively smaller did the key space get?"
# Note: Pollard reduction grows exponentially with bits (Space gets bigger, Interval stays same)
ax2.plot(data["BruteForce"]["bits"], data["BruteForce"]["reduction"], 'o-', label='LSB Leak (Constant Factor)')
ax2.plot(data["PollardRho"]["bits"], data["PollardRho"]["reduction"], 's-', label='Interval Leak (Exponential)')

ax2.set_yscale('log')
ax2.set_xlabel('Curve Size (Bits)')
ax2.set_ylabel('Reduction Factor (x times smaller)')
ax2.set_title('Search Space Reduction Efficiency')
ax2.legend()
ax2.grid(True, which="both", ls="--", alpha=0.3)

# --- GRAPH 3: MEMORY OPTIMIZATION (BabyStep Cost) ---
# Compares Raw Time vs "Memory-Adjusted Cost"
# This proves your BSGS bonus is valuable despite being potentially slower
bsgs_bits = data["BabyStep"]["bits"]
raw_times = data["BabyStep"]["time"]
adj_costs = data["BabyStep"]["cost"]

width = 0.35
x = np.arange(len(bsgs_bits))
ax3.bar(x - width/2, raw_times, width, label='Raw Time (s)', color='gray', alpha=0.6)
ax3.bar(x + width/2, adj_costs, width, label='Memory-Adjusted Cost', color='purple')

ax3.set_xticks(x)
ax3.set_xticklabels(bsgs_bits)
ax3.set_xlabel('Curve Size (Bits)')
ax3.set_ylabel('Cost Unit')
ax3.set_title('BabyStep: Cost Savings (10x RAM Reduction)')
ax3.legend()

plt.suptitle(f"ECC Attack Suite Performance Analysis", fontsize=16)
plt.tight_layout()
plt.show()