# ECDLP Solvers

Implementations of 5 algorithms for solving the Elliptic Curve Discrete Logarithm Problem.

## Quick Start
```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python run_comparisons.py # runs from 10 to 30 bit length
```

## Algorithms

1. **Brute Force** - O(n)
2. **Baby-Step Giant-Step** - O(√n)
3. **Pohlig-Hellman** - O(∑√qᵢ) for smooth orders
4. **Pollard Rho** - O(√n) probabilistic
5. **Las Vegas** - Polynomial probabilistic

## Quick Start

```bash
# Run individual algorithm
python3 <algoname>/main_optimized.py <testcase_path>

# Compare all algorithms
python3 run_comparisons.py <start_bit_size> <end_bit_size> # default : 10 18
# It also generates graphs in graphs/ directory
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

Side-channel attack simulations showing how leaked information breaks ECDLP:

```bash
python3 <algoname>/bonus.py <testcase_path>

python3 run_bonus_comparison.py 
python3 run_bonus_detailed.py
```

**Bonus Results:**
- BruteForce: 16-bit LSB leak → 65,536x search space reduction
- BSGS: 16-bit MSB leak → Reduces √n parameter dramatically
- Pohlig-Hellman: Residue leaks eliminate CRT subproblems

## NOTE :

In run_comparisons.py, uncomment : 
```python
ALGORITHMS = ['BruteForce', 'BabyStep', 'PohligHellman','PollardRho', 'LasVegas']
```