# Smart's Attack for Anomalous Curves

## Overview

Smart's Attack (also known as the Semaev-Smart-Satoh Attack) is a specialized algorithm for solving the Elliptic Curve Discrete Logarithm Problem (ECDLP) on **anomalous curves**.

### What are Anomalous Curves?

A curve E(F_p) is **anomalous** if:
```
#E(F_p) = p
```

Where:
- `#E(F_p)` is the number of points on the curve (order)
- `p` is the prime modulus of the field

### Attack Complexity

- **Time:** O(log p) - effectively instant
- **Space:** O(1)
- **Condition:** Only works when n = p (anomalous)

## Usage

### Basic Usage

```bash
python3 Smart/main.py test_cases/30bit/case_1.txt
```

### With Verification

```bash
python3 Smart/main.py --verify test_cases/35bit/case_2.txt
```

### Optimized Version

```bash
python3 Smart/main_optimized.py --quiet test_cases/40bit/case_3.txt
```

### Run All Test Cases

```bash
oython3 Smart/run.py --start 20 --end 40 --step 5 # optional start end step
```

## Input Format

Test case file format:
```
p
a b
Gx Gy
n
Qx Qy
```

Where:
- `p`: Prime modulus
- `a, b`: Curve parameters (y² = x³ + ax + b)
- `Gx, Gy`: Base point coordinates
- `n`: Order of base point
- `Qx, Qy`: Target point coordinates

## Example Output

```
======================================================================
Smart's Attack - Anomalous Curve Solver
======================================================================
Curve: y² ≡ x³ + 35x + 1 (mod 1073741827)
Base Point G: (123456789, 987654321)
Target Point Q: (456789123, 321987654)
Order n: 1073741827
Field size p: 1073741827
----------------------------------------------------------------------
✓ VULNERABILITY DETECTED: Anomalous Curve
  Condition met: #E(F_p) = p = 1073741827
  Attack: Semaev-Smart-Satoh (p-adic logarithm)
  Complexity: O(log p) - effectively instant
----------------------------------------------------------------------
Executing attack...

✓ ATTACK SUCCESSFUL
Solution: d = 123456
Time: 0.000234s
Verification: 123456*G = Q ✓
======================================================================
```

## Files

- `main.py` - Full implementation with detailed output
- `main_optimized.py` - Optimized version with minimal output
- `smart.py` - Library module for use in other scripts
- `gen.py` - Test case generator for anomalous curves

## Algorithm Details

### Mathematical Background

On an anomalous curve, the ECDLP can be solved using:

1. **Hensel Lifting:** Lift points from F_p to the p-adic field Q_p
2. **p-adic Logarithm:** Compute discrete logs in the formal group
3. **Modular Division:** Solve for d in Z/pZ

### Implementation Notes

This implementation provides:
- ✓ Anomalous condition detection
- ✓ Input validation and error handling
- ✓ Solution verification
- ✓ Progress indicators for large searches
- ✓ Answer file verification (--verify flag)

For production use, integrate p-adic arithmetic libraries like:
- SageMath's p-adic module
- PARI/GP for computational number theory

## Security Implications

**Never use anomalous curves in production!**

If a curve is anomalous:
- All private keys can be recovered instantly
- No computational security remains
- The ECDLP is trivially solvable

Standard curves (NIST P-256, secp256k1, etc.) are specifically chosen to **not** be anomalous.

## References

1. Nigel Smart (1999): "The Discrete Logarithm Problem on Elliptic Curves of Trace One"
2. Igor Semaev (1998): "Evaluation of discrete logarithms in a group of p-torsion points"
3. Takakazu Satoh (2000): "The canonical lift of an ordinary elliptic curve"

## Testing

Run all test cases:
```bash
for case in test_cases/30bit/case_*.txt; do
    echo "Testing $case..."
    python3 Smart/main.py "$case"
done
```
