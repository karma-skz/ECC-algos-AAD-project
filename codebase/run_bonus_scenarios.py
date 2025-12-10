#!/usr/bin/env python3
"""
MASTER BONUS RUNNER
Runs all bonus scenarios on a specific test case and generates a summary.
Usage: python3 run_bonus_scenarios.py <path_to_test_case>
"""
import subprocess
import sys
import json
import time
from pathlib import Path
from datetime import datetime

def print_separator(char="=", length=80):
    print(char * length)

def run_scenario(algo, test_file):
    script = Path(algo) / 'bonus.py'
    if not script.exists(): 
        print(f"⚠ Script not found: {script}")
        return None
    
    print(f"\n>>> Running {algo} Bonus...")
    try:
        # Capture stdout to parse the result
        result = subprocess.run(
            [sys.executable, str(script), str(test_file)], 
            capture_output=True, 
            text=True,
            check=False
        )
        
        # Print the output to console so user sees progress
        print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)

        # Parse the last line that starts with BONUS_RESULT
        for line in result.stdout.splitlines():
            if line.startswith("BONUS_RESULT:"):
                try:
                    json_str = line.replace("BONUS_RESULT:", "").strip()
                    return json.loads(json_str)
                except json.JSONDecodeError:
                    print("⚠ Failed to parse bonus result JSON")
                    return None
        
        return None

    except KeyboardInterrupt:
        print("\n[Aborted by user]")
        sys.exit(1)
    except Exception as e:
        print(f"Error running {algo}: {e}")
        return None

def generate_html_report(results, case_path, bits):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Prepare data for charts
    algos = [r['algo'] for r in results]
    # Ensure non-zero values for chart visibility (min 0.000001)
    times = [max(r['time'], 0.000001) for r in results]
    statuses = [r['status'] for r in results]
    colors = ['#4CAF50' if s == 'success' else '#F44336' for s in statuses]
    
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>ECC Bonus Scenarios Report</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 20px; background-color: #f5f5f5; }}
        .container {{ max_width: 1000px; margin: 0 auto; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }}
        h1 {{ color: #333; border-bottom: 2px solid #eee; padding-bottom: 10px; }}
        .meta {{ color: #666; margin-bottom: 20px; }}
        .card {{ background: #fff; border: 1px solid #ddd; border-radius: 4px; padding: 15px; margin-bottom: 20px; }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 20px; }}
        th, td {{ padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background-color: #f8f9fa; }}
        .success {{ color: green; font-weight: bold; }}
        .fail {{ color: red; font-weight: bold; }}
        .chart-container {{ position: relative; height: 400px; width: 100%; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>ECC Side-Channel & Bonus Scenarios Report</h1>
        <div class="meta">
            <p><strong>Date:</strong> {timestamp}</p>
            <p><strong>Test Case:</strong> {case_path}</p>
            <p><strong>Curve Size:</strong> {bits}-bit</p>
        </div>

        <div class="card">
            <h2>Performance Overview</h2>
            <div class="chart-container">
                <canvas id="perfChart"></canvas>
            </div>
        </div>

        <div class="card">
            <h2>Detailed Results</h2>
            <table>
                <thead>
                    <tr>
                        <th>Algorithm</th>
                        <th>Status</th>
                        <th>Time (s)</th>
                        <th>Steps/Ops</th>
                        <th>Details</th>
                    </tr>
                </thead>
                <tbody>
    """
    
    for r in results:
        status_class = "success" if r['status'] == 'success' else "fail"
        details_str = ", ".join([f"{k}={v}" for k,v in r['details'].items()])
        html_content += f"""
                    <tr>
                        <td>{r['algo']}</td>
                        <td class="{status_class}">{r['status'].upper()}</td>
                        <td>{r['time']:.6f}</td>
                        <td>{r['steps']:,}</td>
                        <td>{details_str}</td>
                    </tr>
        """

    html_content += f"""
                </tbody>
            </table>
        </div>
    </div>

    <script>
        const ctx = document.getElementById('perfChart').getContext('2d');
        new Chart(ctx, {{
            type: 'bar',
            data: {{
                labels: {json.dumps(algos)},
                datasets: [{{
                    label: 'Execution Time (seconds)',
                    data: {json.dumps(times)},
                    backgroundColor: {json.dumps(colors)},
                    borderWidth: 1
                }}]
            }},
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                scales: {{
                    y: {{
                        beginAtZero: true,
                        title: {{
                            display: true,
                            text: 'Time (s)'
                        }}
                    }}
                }},
                plugins: {{
                    legend: {{
                        display: false
                    }},
                    title: {{
                        display: true,
                        text: 'Execution Time by Scenario (Lower is Better)'
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
    output_file = output_dir / "report.html"
    
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(html_content)
    
    print(f"\n[+] Report generated: {output_file.absolute()}")
    return output_file

def main():
    # 1. Parse Arguments
    if len(sys.argv) > 1:
        case_path = Path(sys.argv[1])
    else:
        # Default fallback
        print("No test case specified. Looking for a default...")
        case_path = Path("test_cases/20bit/case_1.txt")
        if not case_path.exists():
            # Find first available txt file
            cases = list(Path("test_cases").glob("*/*.txt"))
            if cases:
                case_path = cases[0]
            else:
                print("Error: No test cases found in test_cases/ folder.")
                sys.exit(1)

    if not case_path.exists():
        print(f"Error: File not found: {case_path}")
        sys.exit(1)

    # 2. Print Header
    print_separator("=")
    print(f"  ECC SIDE-CHANNEL / BONUS SCENARIOS TEST SUITE")
    print_separator("=")
    print(f"Target Case: {case_path}")
    
    # Check bit size from path name (heuristic)
    bits = 0
    try:
        bits = int(case_path.parent.name.replace('bit', ''))
        print(f"Bit Length:  {bits}-bit curve")
    except:
        pass
    
    print_separator("=")

    # 3. Run All Scenarios
    algos = ['BruteForce', 'BabyStep', 'PollardRho', 'PohligHellman', 'LasVegas']
    results = []
    
    for algo in algos:
        res = run_scenario(algo, case_path)
        if res:
            results.append(res)
        else:
            # Add a dummy failed result if script failed or didn't output JSON
            results.append({
                "algo": algo,
                "status": "error",
                "time": 0,
                "steps": 0,
                "details": {"error": "Script execution failed"}
            })

    print("\n")
    print_separator("=")
    print("[+] SCENARIO SUITE COMPLETE")
    print_separator("=")
    
    # 4. Generate Report
    if results:
        report_path = generate_html_report(results, case_path, bits)
        # Try to open the report
        try:
            if sys.platform == 'win32':
                subprocess.run(['start', str(report_path)], shell=True)
            elif sys.platform == 'darwin':
                subprocess.run(['open', str(report_path)])
            else:
                subprocess.run(['xdg-open', str(report_path)])
        except:
            pass

if __name__ == "__main__":
    main()