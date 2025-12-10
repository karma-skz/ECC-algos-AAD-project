#!/usr/bin/env python3
"""
AGGREGATE BONUS RUNNER
Runs bonus scenarios across ALL test cases to generate a performance curve.
Usage: python3 run_aggregate_bonus.py
"""
import sys
import json
import time
from pathlib import Path
from datetime import datetime
import run_bonus_scenarios  # Import the single-runner logic

def get_test_cases():
    """Find one test case per bit size, sorted by size."""
    cases = []
    base_dir = Path("test_cases")
    
    if not base_dir.exists():
        print("Error: test_cases directory not found.")
        return []

    for d in base_dir.iterdir():
        if d.is_dir() and d.name.endswith("bit"):
            try:
                bits = int(d.name.replace("bit", ""))
                # Prefer case_1.txt, but take any txt if missing
                case_file = d / "case_1.txt"
                if not case_file.exists():
                    txt_files = list(d.glob("*.txt"))
                    if txt_files:
                        case_file = txt_files[0]
                    else:
                        continue
                
                cases.append({
                    "bits": bits,
                    "path": case_file
                })
            except ValueError:
                continue
    
    # Sort by bit size
    return sorted(cases, key=lambda x: x["bits"])

def generate_aggregate_report(all_results):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Extract data for charts
    # Structure: { "BruteForce": { 10: 0.001, 11: 0.002 ... }, ... }
    chart_data = {}
    bit_sizes = sorted(list({r["bits"] for r in all_results}))
    
    for r in all_results:
        algo = r["algo"]
        bits = r["bits"]
        time_val = max(r["time"], 0.000001) # Ensure visibility
        
        if algo not in chart_data:
            chart_data[algo] = []
        
        # We need to align data with bit_sizes
        # This simple append assumes we process in order, which we do.
        chart_data[algo].append(time_val)

    # Colors for lines
    colors = {
        "BruteForce": "#FF6384",
        "BabyStep": "#36A2EB",
        "PollardRho": "#FFCE56",
        "PohligHellman": "#4BC0C0",
        "LasVegas": "#9966FF"
    }

    datasets = []
    for algo, data in chart_data.items():
        datasets.append({
            "label": algo,
            "data": data,
            "borderColor": colors.get(algo, "#000000"),
            "fill": False,
            "tension": 0.1
        })

    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>ECC Bonus Scenarios - Aggregate Performance</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 20px; background-color: #f5f5f5; }}
        .container {{ max_width: 1200px; margin: 0 auto; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }}
        h1 {{ color: #333; border-bottom: 2px solid #eee; padding-bottom: 10px; }}
        .meta {{ color: #666; margin-bottom: 20px; }}
        .card {{ background: #fff; border: 1px solid #ddd; border-radius: 4px; padding: 15px; margin-bottom: 20px; }}
        .chart-container {{ position: relative; height: 600px; width: 100%; }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 20px; font-size: 0.9em; }}
        th, td {{ padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background-color: #f8f9fa; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>ECC Bonus Scenarios - Aggregate Performance Curve</h1>
        <div class="meta">
            <p><strong>Date:</strong> {timestamp}</p>
            <p><strong>Range:</strong> {bit_sizes[0]}-bit to {bit_sizes[-1]}-bit</p>
        </div>

        <div class="card">
            <h2>Execution Time vs Curve Size</h2>
            <p>Note: Y-axis is logarithmic to handle wide variations in timing.</p>
            <div class="chart-container">
                <canvas id="aggChart"></canvas>
            </div>
        </div>

        <div class="card">
            <h2>Raw Data</h2>
            <table>
                <thead>
                    <tr>
                        <th>Bits</th>
                        <th>Algorithm</th>
                        <th>Time (s)</th>
                        <th>Status</th>
                    </tr>
                </thead>
                <tbody>
    """
    
    for r in all_results:
        html_content += f"""
                    <tr>
                        <td>{r['bits']}</td>
                        <td>{r['algo']}</td>
                        <td>{r['time']:.6f}</td>
                        <td>{r['status']}</td>
                    </tr>
        """

    html_content += f"""
                </tbody>
            </table>
        </div>
    </div>

    <script>
        const ctx = document.getElementById('aggChart').getContext('2d');
        new Chart(ctx, {{
            type: 'line',
            data: {{
                labels: {json.dumps(bit_sizes)},
                datasets: {json.dumps(datasets)}
            }},
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                scales: {{
                    y: {{
                        type: 'logarithmic',
                        title: {{
                            display: true,
                            text: 'Time (seconds) - Log Scale'
                        }},
                        ticks: {{
                            callback: function(value, index, values) {{
                                return Number(value.toString());
                            }}
                        }}
                    }},
                    x: {{
                        title: {{
                            display: true,
                            text: 'Curve Size (Bits)'
                        }}
                    }}
                }},
                plugins: {{
                    title: {{
                        display: true,
                        text: 'Algorithm Performance Scaling'
                    }},
                    tooltip: {{
                        callbacks: {{
                            label: function(context) {{
                                return context.dataset.label + ': ' + context.raw.toFixed(6) + ' s';
                            }}
                        }}
                    }}
                }}
            }}
        }});
    </script>
</body>
</html>
    """
    
    output_dir = Path("bonus_results")
    output_dir.mkdir(exist_ok=True)
    output_file = output_dir / "aggregate_report.html"
    
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(html_content)
    
    print(f"\n[+] Aggregate Report generated: {output_file.absolute()}")
    return output_file

def main():
    print(">>> Finding test cases...")
    cases = get_test_cases()
    if not cases:
        print("No test cases found.")
        return

    print(f"Found {len(cases)} test cases: {[c['bits'] for c in cases]}")
    
    all_results = []
    algos = ['BruteForce', 'BabyStep', 'PollardRho', 'PohligHellman', 'LasVegas']

    for c in cases:
        bits = c['bits']
        path = c['path']
        print(f"\n{'#'*60}")
        print(f"PROCESSING {bits}-BIT CURVE")
        print(f"{'#'*60}")
        
        for algo in algos:
            # Skip very large curves for BabyStep if it's going to be too slow
            # BabyStep bonus is O(sqrt(N/100)). For 40 bits, sqrt(2^40/100) ~ 100k. Fast.
            # For 50 bits, sqrt(2^50/100) ~ 3M. Might take a few seconds. Acceptable.
            
            res = run_bonus_scenarios.run_scenario(algo, path)
            if res:
                res['bits'] = bits
                all_results.append(res)
            else:
                all_results.append({
                    "algo": algo,
                    "bits": bits,
                    "status": "error",
                    "time": 0,
                    "steps": 0,
                    "details": {}
                })

    print("\n" + "="*60)
    print("AGGREGATION COMPLETE")
    print("="*60)
    
    report_path = generate_aggregate_report(all_results)
    
    # Try to open
    try:
        if sys.platform == 'win32':
            import subprocess
            subprocess.run(['start', str(report_path)], shell=True)
        elif sys.platform == 'darwin':
            subprocess.run(['open', str(report_path)]) #type:ignore
        else:
            subprocess.run(['xdg-open', str(report_path)]) #type:ignore
    except:
        pass

if __name__ == "__main__":
    main()
