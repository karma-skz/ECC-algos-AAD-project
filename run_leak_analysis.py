#!/usr/bin/env python3
"""
LEAK ANALYSIS RUNNER (MULTI-BIT)
Runs specific leak scenarios across MULTIPLE bit lengths to generate comparative graphs.
"""
import subprocess
import sys
import json
import time
import argparse
from pathlib import Path
from datetime import datetime
from utils import load_input

# Configuration
DEFAULT_BIT_LENGTHS = [15, 20, 25, 30, 35, 40]
DEFAULT_LEAKS = [0, 4, 8, 12, 16, 20, 24]
TIMEOUT_SEC = 10  # Fast timeout to skip impossible cases

def get_test_case(bits, case_num):
    """Finds a valid test case for the requested bit length and case number."""
    # Try exact folder match first
    folder = Path(f"test_cases/{bits}bit")
    
    # Helper to check common patterns
    patterns = [
        f"case_{case_num}.txt",
        f"test_{case_num}.txt",
        f"testcase_{case_num}.txt",
        f"*{case_num}.txt" # Catch-all
    ]

    if folder.exists():
        for pat in patterns:
            matches = list(folder.glob(pat))
            if matches: return matches[0]
    
    # Fallback: search all folders for bit string and case number
    for p in Path("test_cases").glob(f"*/*{case_num}.txt"):
        if f"{bits}bit" in str(p):
            return p
            
    return None

def run_algo_scenario(algo, test_file, args):
    script = Path(algo) / 'bonus.py'
    if not script.exists(): return None
    
    cmd = [sys.executable, str(script), str(test_file)] + args
    try:
        # Use a timeout to prevent hanging on hard cases
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_SEC)
        for line in result.stdout.splitlines():
            if line.startswith("BONUS_RESULT:"):
                return json.loads(line.replace("BONUS_RESULT:", "").strip())
    except subprocess.TimeoutExpired:
        return {"algo": algo, "status": "timeout", "time": float(TIMEOUT_SEC), "details": {}}
    except Exception as e:
        print(f"Error running {algo}: {e}")
    return None

def run_suite(target_bits, case_num, target_leaks):
    # Data structure: data[algo][bit_len] = [ {x: leak, y: time}, ... ]
    suite_data = {
        "BruteForce": {},
        "BabyStep": {},
        "PollardRho": {},
        "Kangaroo": {},
        "LasVegas": {}
    }

    print(f"\n{'='*60}")
    print(f"STARTING MULTI-BIT LEAK ANALYSIS")
    print(f"Target Sizes: {target_bits}")
    print(f"Target Leaks: {target_leaks}")
    print(f"Test Case: #{case_num}")
    print(f"{'='*60}")

    for bits in target_bits:
        case_path = get_test_case(bits, case_num)
        if not case_path:
            print(f"\n[!] Skipping {bits}-bit: No test case #{case_num} found.")
            continue
            
        try:
            _, _, _, _, n, _ = load_input(case_path)
            real_bit_len = n.bit_length()
            print(f"\n>>> Testing {bits}-bit folder (File: {case_path.name}, N is {real_bit_len}-bit)...")
        except:
            print(f"\n[!] Error loading {case_path}")
            continue

        # --- 1. LSB Leakage (BruteForce, BSGS, Rho) ---
        # Filter leaks that are larger than N
        leaks = [l for l in target_leaks if l < real_bit_len]

        for algo in ["BruteForce", "BabyStep", "PollardRho"]:
            if bits not in suite_data[algo]: suite_data[algo][bits] = []
            
            for leak in leaks:
                print(f"  [{algo}] Leak {leak} bits...", end="", flush=True)
                res = run_algo_scenario(algo, case_path, ["--leak-bits", str(leak)])
                
                if res:
                    t = res['time']
                    status = res['status']
                    print(f" {t:.4f}s ({status})")
                    
                    if status == 'success':
                        suite_data[algo][bits].append({"x": leak, "y": t})
                    elif status == 'timeout':
                        suite_data[algo][bits].append({"x": leak, "y": TIMEOUT_SEC})
                else:
                    print(" Error")

        # --- 2. Interval (Kangaroo) ---
        if bits not in suite_data["Kangaroo"]: suite_data["Kangaroo"][bits] = []
        
        widths = []
        curr = n
        while curr > 100:
            widths.append(curr)
            curr //= 16 # Faster reduction for speed
            if len(widths) > 5: break
            
        for w in widths:
            print(f"  [Kangaroo] Width {w}...", end="", flush=True)
            res = run_algo_scenario("PollardRho", case_path, ["--interval-width", str(w)])
            if res:
                print(f" {res['time']:.4f}s")
                if res['status'] == 'success':
                    factor = n / w
                    suite_data["Kangaroo"][bits].append({"x": factor, "y": res['time']})
            else: print(" Error")

        # --- 3. Approx (Las Vegas) ---
        if bits not in suite_data["LasVegas"]: suite_data["LasVegas"][bits] = []
        
        errors = [100, 1000, 10000, 50000]
        for err in errors:
            print(f"  [LasVegas] Error {err}...", end="", flush=True)
            res = run_algo_scenario("LasVegas", case_path, ["--approx-error", str(err)])
            if res:
                print(f" {res['time']:.4f}s")
                if res['status'] == 'success':
                    suite_data["LasVegas"][bits].append({"x": err, "y": res['time']})
            else: print(" Error")

    return suite_data

def generate_html(data, target_bits):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Helper to generate datasets for a chart
    def make_datasets(algo_name):
        datasets = []
        colors = ['#FF6384', '#36A2EB', '#FFCE56', '#4BC0C0', '#9966FF', '#FF9F40', '#E7E9ED', '#71B37C']
        
        for i, bits in enumerate(target_bits):
            if bits not in data[algo_name] or not data[algo_name][bits]:
                continue
                
            points = data[algo_name][bits]
            # Sort by X
            points.sort(key=lambda p: p['x'])
            
            datasets.append({
                "label": f"{bits}-bit Curve",
                "data": points,
                "borderColor": colors[i % len(colors)],
                "fill": False,
                "tension": 0.1
            })
        return datasets

    # Pre-calculate JSON strings for datasets
    bf_data_json = json.dumps(make_datasets("BruteForce"))
    bsgs_data_json = json.dumps(make_datasets("BabyStep"))
    rho_data_json = json.dumps(make_datasets("PollardRho"))
    kang_data_json = json.dumps(make_datasets("Kangaroo"))
    lv_data_json = json.dumps(make_datasets("LasVegas"))

    html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>ECC Multi-Bit Leak Analysis</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body {{ font-family: 'Segoe UI', sans-serif; margin: 20px; background: #f5f7fa; color: #333; }}
        .container {{ max_width: 1200px; margin: 0 auto; }}
        .header {{ background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 20px; }}
        .grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
        .chart-card {{ background: white; padding: 15px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        h1 {{ margin: 0 0 10px 0; font-size: 24px; }}
        h2 {{ font-size: 18px; color: #555; border-bottom: 1px solid #eee; padding-bottom: 10px; }}
        .full-width {{ grid-column: 1 / -1; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>ECC Side-Channel Analysis: Multi-Bit Comparison</h1>
            <p>Generated: {timestamp} | Timeout: {TIMEOUT_SEC}s</p>
            <p>This report compares how different algorithms perform under partial information leakage across increasing curve sizes.</p>
        </div>

        <div class="grid">
            <!-- 1. Brute Force -->
            <div class="chart-card">
                <h2>1. Brute Force: Impact of LSB Leakage</h2>
                <canvas id="bfChart"></canvas>
                <p><small>X: Bits Leaked | Y: Time (s)</small></p>
            </div>

            <!-- 2. Baby Step -->
            <div class="chart-card">
                <h2>2. Baby Step Giant Step: Impact of LSB Leakage</h2>
                <canvas id="bsgsChart"></canvas>
                <p><small>X: Bits Leaked | Y: Time (s)</small></p>
            </div>

            <!-- 3. Pollard Rho -->
            <div class="chart-card">
                <h2>3. Pollard Rho: Impact of LSB Leakage</h2>
                <canvas id="rhoChart"></canvas>
                <p><small>X: Bits Leaked | Y: Time (s)</small></p>
            </div>

            <!-- 4. Kangaroo -->
            <div class="chart-card">
                <h2>4. Pollard's Kangaroo: Interval Reduction</h2>
                <canvas id="kangChart"></canvas>
                <p><small>X: Reduction Factor (N/Width) | Y: Time (s)</small></p>
            </div>

            <!-- 5. Las Vegas -->
            <div class="chart-card full-width">
                <h2>5. Las Vegas: Approximate Key Knowledge</h2>
                <canvas id="lvChart"></canvas>
                <p><small>X: Error Margin (+/-) | Y: Time (s)</small></p>
            </div>
        </div>
    </div>

    <script>
        const commonOptions = {{
            scales: {{
                y: {{ type: 'logarithmic', title: {{ display: true, text: 'Time (s)' }} }},
                x: {{ title: {{ display: true, text: 'Parameter' }} }}
            }},
            plugins: {{ legend: {{ position: 'right' }} }}
        }};

        // 1. Brute Force
        new Chart(document.getElementById('bfChart'), {{
            type: 'line',
            data: {{ datasets: {bf_data_json} }},
            options: {{ 
                ...commonOptions, 
                scales: {{ 
                    ...commonOptions.scales, 
                    x: {{ 
                        type: 'linear',
                        position: 'bottom',
                        title: {{ display: true, text: 'Bits Leaked' }} 
                    }} 
                }} 
            }}
        }});

        // 2. BSGS
        new Chart(document.getElementById('bsgsChart'), {{
            type: 'line',
            data: {{ datasets: {bsgs_data_json} }},
            options: {{ 
                ...commonOptions, 
                scales: {{ 
                    ...commonOptions.scales, 
                    x: {{ 
                        type: 'linear',
                        position: 'bottom',
                        title: {{ display: true, text: 'Bits Leaked' }} 
                    }} 
                }} 
            }}
        }});

        // 3. Rho
        new Chart(document.getElementById('rhoChart'), {{
            type: 'line',
            data: {{ datasets: {rho_data_json} }},
            options: {{ 
                ...commonOptions, 
                scales: {{ 
                    ...commonOptions.scales, 
                    x: {{ 
                        type: 'linear',
                        position: 'bottom',
                        title: {{ display: true, text: 'Bits Leaked' }} 
                    }} 
                }} 
            }}
        }});

        // 4. Kangaroo
        new Chart(document.getElementById('kangChart'), {{
            type: 'line',
            data: {{ datasets: {kang_data_json} }},
            options: {{ 
                ...commonOptions, 
                scales: {{ 
                    ...commonOptions.scales, 
                    x: {{ 
                        type: 'logarithmic', 
                        title: {{ display: true, text: 'Reduction Factor (N/w)' }} 
                    }} 
                }} 
            }}
        }});

        // 5. Las Vegas
        new Chart(document.getElementById('lvChart'), {{
            type: 'line',
            data: {{ datasets: {lv_data_json} }},
            options: {{ 
                ...commonOptions, 
                scales: {{ 
                    ...commonOptions.scales, 
                    x: {{ 
                        type: 'linear', 
                        title: {{ display: true, text: 'Error Margin' }} 
                    }} 
                }} 
            }}
        }});
    </script>
</body>
</html>
    """
    
    # Create the directory if it doesn't exist
    output_dir = Path("bonus_results")
    output_dir.mkdir(exist_ok=True)

    # Update path to save inside that directory
    out_path = output_dir / "multi_bit_leak_report.html"
    
    with open(out_path, "w") as f:
        f.write(html)
    print(f"\n[+] Report generated: {out_path.absolute()}")
    return out_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run multi-bit leak analysis.")
    parser.add_argument("--bits", type=int, nargs="+", default=DEFAULT_BIT_LENGTHS, help="List of bit lengths to test (e.g. 15 20 30)")
    parser.add_argument("--case", type=int, default=1, help="Test case number to use (1-5)")
    parser.add_argument("--leaks", type=int, nargs="+", default=DEFAULT_LEAKS, help="List of leak bits to test (e.g. 0 4 8 12)")
    args = parser.parse_args()

    data = run_suite(args.bits, args.case, args.leaks)
    report = generate_html(data, args.bits)
    if sys.platform == 'darwin':
        subprocess.run(['open', str(report)])
