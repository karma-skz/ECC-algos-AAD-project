# ECC Algorithms for ECDLP - AAD Project

A modular, well-structured implementation of various algorithms for solving the Elliptic Curve Discrete Logarithm Problem (ECDLP).

## Project Structure

```
.
├── utils/                      # Shared utility modules
│   ├── __init__.py
│   ├── ecc_utils.py           # EllipticCurve class and point operations
│   ├── mod_utils.py           # Modular arithmetic utilities
│   ├── io_utils.py            # Input/output handling
│   └── generator.py           # Test case generator
├── BruteForce/                # Brute force algorithm
│   ├── main.py
│   └── input/
├── BabyStep/                  # Baby-Step Giant-Step (BSGS)
│   ├── main.py
│   └── input/
├── PohligHellman/             # Pohlig-Hellman algorithm
│   ├── main.py
│   └── input/
├── PollardRho/                # Pollard's Rho algorithm
│   ├── main.py
│   └── input/
└── LasVegas/                  # Las Vegas probabilistic algorithm
    ├── main.py
    └── input/
```

## Algorithms Implemented

### 1. Brute Force
- **Complexity:** O(n) time, O(1) space
- **Method:** Exhaustive search through all possible values
- **Best for:** Very small orders (n < 10,000)

### 2. Baby-Step Giant-Step (BSGS)
- **Complexity:** O(√n) time, O(√n) space
- **Method:** Meet-in-the-middle approach using hash table
- **Best for:** Small to medium orders (n < 10^12)

### 3. Pohlig-Hellman
- **Complexity:** O(∑ eᵢ(log n + √qᵢ)) where n = q₁^e₁ × ... × qₖ^eₖ
- **Method:** Exploits smooth order (product of small primes)
- **Best for:** Cases where n has small prime factors

### 4. Pollard's Rho
- **Complexity:** O(√n) expected time, O(1) space
- **Method:** Probabilistic collision detection using Floyd's cycle
- **Best for:** General purpose, space-efficient

### 5. Las Vegas
- **Complexity:** Polynomial per attempt (probabilistic)
- **Method:** Linear algebra over finite fields with summation polynomials
- **Best for:** Experimental/research purposes

## Shared Utilities

### EllipticCurve Class (`utils/ecc_utils.py`)
```python
from utils import EllipticCurve

curve = EllipticCurve(a=497, b=1768, p=9739)
G = (1804, 5368)
Q = (3138, 1774)

# Point operations
R = curve.add(G, Q)                    # Point addition
S = curve.scalar_multiply(5, G)       # Scalar multiplication
T = curve.negate(G)                    # Point negation
valid = curve.is_on_curve(G)          # Validation
```

### Modular Arithmetic (`utils/mod_utils.py`)
```python
from utils import mod_inv, extended_gcd, crt_combine

# Modular inverse
inv = mod_inv(3, 7)  # 3^(-1) mod 7 = 5

# Extended GCD
gcd, x, y = extended_gcd(12, 18)  # gcd=6, 12x + 18y = 6

# Chinese Remainder Theorem
congruences = [(2, 3), (3, 5), (2, 7)]  # x ≡ 2 (mod 3), x ≡ 3 (mod 5), x ≡ 2 (mod 7)
result, modulus = crt_combine(congruences)
```

### Input/Output (`utils/io_utils.py`)
```python
from utils import load_input

# Load ECDLP problem from file
p, a, b, G, n, Q = load_input(Path("input/test_1.txt"))
```

## Input File Format

All algorithms use the same 5-line input format:

```
<p>              # Prime modulus
<a> <b>          # Curve coefficients (y² = x³ + ax + b)
<Gx> <Gy>        # Base point G coordinates
<n>              # Order of base point
<Qx> <Qy>        # Target point Q coordinates
```

**Example:**
```
9739
497 1768
1804 5368
9735
3138 1774
```

This represents:
- Curve: y² = x³ + 497x + 1768 (mod 9739)
- Base point: G = (1804, 5368)
- Order: n = 9735
- Target: Q = (3138, 1774)
- Find: d such that Q = d·G

## Usage

### Running Individual Algorithms

```bash
# Brute Force
cd BruteForce
python3 main.py input/test_1.txt

# Baby-Step Giant-Step
cd BabyStep
python3 main.py input/test_1.txt

# Pohlig-Hellman
cd PohligHellman
python3 main.py input/test_1.txt

# Pollard Rho
cd PollardRho
python3 main.py input/test_1.txt

# Las Vegas
cd LasVegas
python3 main.py input/test_1.txt
```

### Generating Test Cases

```bash
cd utils
python3 generator.py
```

This generates 5 test cases in the `input/` directory with corresponding answer files.

## Key Features

### ✅ Modular Design
- Shared utilities eliminate code duplication
- Clean separation of concerns
- Easy to maintain and extend

### ✅ Professional Code Quality
- Comprehensive docstrings
- Type hints throughout
- Clear variable names
- Proper error handling

### ✅ Verified Correctness
- All solutions are verified: Q = d·G
- Extensive validation of inputs
- Clear success/failure reporting

### ✅ Performance Optimized
- Efficient double-and-add for scalar multiplication
- Optimized RREF for linear algebra
- Memory-efficient data structures

### ✅ Educational Value
- Clear algorithm descriptions in comments
- Time/space complexity documented
- Easy to understand implementations

## Algorithm Comparison

| Algorithm | Time | Space | Best Use Case |
|-----------|------|-------|---------------|
| Brute Force | O(n) | O(1) | n < 10⁴ |
| BSGS | O(√n) | O(√n) | n < 10¹² |
| Pohlig-Hellman | O(∑√qᵢ) | O(√max(qᵢ)) | Smooth n |
| Pollard Rho | O(√n) | O(1) | General |
| Las Vegas | O(poly) | O(n'²) | Research |

## Dependencies

- Python 3.7+
- No external libraries required (uses only standard library)

## Testing

All algorithms have been tested with the provided test cases:

```bash
# Test Brute Force
cd BruteForce && python3 main.py input/test_1.txt

# Test BSGS  
cd BabyStep && python3 main.py input/test_1.txt

# Test Pohlig-Hellman
cd PohligHellman && python3 main.py input/test_1.txt

# Test Pollard Rho
cd PollardRho && python3 main.py input/test_1.txt

# Test Las Vegas
cd LasVegas && python3 main.py input/test_1.txt
```

**Expected Output:**
```
Solving ECDLP using [Algorithm Name]...
Curve: y^2 = x^3 + [a]x + [b] (mod [p])
G = ([Gx], [Gy]), Q = ([Qx], [Qy]), n = [n]

==================================================
Solution: d = [value]
Time: [seconds] seconds
Verification: PASSED
==================================================
```

## Implementation Details

### Point Representation
Points are represented as tuples `(x, y)`, with `None` representing the point at infinity (identity element).

### Scalar Multiplication
Uses the efficient double-and-add algorithm:
```python
def scalar_multiply(k, P):
    R = O  # Identity
    while k > 0:
        if k & 1: R = R + P
        P = 2P
        k >>= 1
    return R
```

### Error Handling
- Validates curve is non-singular (discriminant ≠ 0)
- Checks points lie on the curve
- Verifies order condition n·G = O
- Handles edge cases (point at infinity, etc.)

## Academic Context

This project implements algorithms discussed in:
- Koblitz, "Elliptic Curve Cryptosystems"
- Silverman & Tate, "Rational Points on Elliptic Curves"
- Shanks, "Class Number, a Theory of Factorization and Genera"
- Pollard, "Monte Carlo Methods for Index Computation"

## License

Academic project for AAD course, Semester 2.

## Authors

- Course: Advanced Algorithm Design (AAD)
- Semester: 2-1
- Institution: [Your Institution]

---

**Note:** These implementations are for educational purposes. For production cryptographic applications, use established libraries like `cryptography` or `PyCryptodome`.
