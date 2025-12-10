# ECDLP Solvers

Check ``FINAL_SUBMISSIONS`` folder for project report and final presentation.
This is the codebase we developed for the same.

Implementations of 5 algorithms for solving the Elliptic Curve Discrete Logarithm Problem.

## Quick Start
```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt

cd codebase
python run_comparison.py 10 30  # runs from 10 to 30 bit length
```

## Algorithms

1. **Brute Force** - O(n)
2. **Baby-Step Giant-Step** - O(√n)
3. **Pohlig-Hellman** - O(∑√qᵢ) for smooth orders
4. **Pollard Rho** - O(√n) probabilistic
5. **Las Vegas** - Polynomial probabilistic

## Quick Start

```bash
cd codebase

python3 generate_test_cases.py <start_bit> <end_bit>  

python3 algo_name/main_optimized.py <testcase_path>  # test specific algo

python3 run_comparison.py <start_bit> <end_bit>  # compare all algos
# make sure to install matplotlib before this for graphs

# Optional demo (if present)
python3 demo.py <testcase_path>  # demo with user input
```

## Input Format

```
p           # Prime modulus
a b         # Curve coefficients y² = x³ + ax + b
Gx Gy       # Base point G
n           # Order of G
Qx Qy       # Target point Q
```

Output: Secret `d` where Q = d·G

## Bonus: Partial Key Leakage

```bash
cd codebase

python3 <algoname>/bonus.py <testcase_path>

python3 run_bonus_scenarios.py <testcase_path>
```

## Project Structure

- codebase/ – all algorithm implementations, scripts, utilities, test cases, graphs, and web interface.
- FINAL_SUBMISSIONS/ – final_report.pdf and presentation.pdf (final submitted artifacts only).
- requirements.txt, README.md – kept at the repository root for environment setup and documentation.

## NOTE :

1. In run_comparisons.py, keep the required algos : 
    ```python
    ALGORITHMS = ['BruteForce', 'BabyStep', 'PohligHellman','PollardRho', 'LasVegas']
    ```
2. For rigorous tests, U may need to uncomment :
    ```python
    # --- SMART LIMITS ---
    if algo == 'BruteForce' and bits > 24:
        print(f"  {algo:15s}: SKIPPED (exponential time)")
        continue
    if algo == 'BabyStep' and bits > 50:
        print(f"  {algo:15s}: SKIPPED (memory limit)")
        continue
    ```
3. In LasVegas/main_optimized.py, u can change the limit :

    ```python
    limit = min(max_attempts, max(2000, expected_points * 10))  # 10 -> other smaller no.
    ```
