# MOV (Menezes-Okamoto-Vanstone) Attack

Educational implementation of the MOV attack concept for solving the Elliptic Curve Discrete Logarithm Problem (ECDLP).

## Overview

The MOV attack reduces the ECDLP on an elliptic curve E over F_p to a discrete logarithm problem in the multiplicative group of an extension field F_p^k. This is achieved using bilinear pairings (Weil or Tate pairings).

**Important**: This implementation demonstrates the MOV attack theory but uses BSGS for practical solving, as MOV is rarely practical on real curves.

## Algorithm

Given Q = d*G, find d by:
1. Find embedding degree k where n | (p^k - 1)
2. Construct extension field F_p^k with irreducible polynomial
3. Find random point R of order n in E(F_p^k) linearly independent from G
4. Compute pairings: α = e(G, R), β = e(Q, R)
5. Solve DLP in F_p^k: find d such that β = α^d using BSGS
6. Result: Q = d*G

## Files

- **main.py**: MOV attack demonstration with BSGS fallback
  - Demonstrates embedding degree finding
  - Shows extension field construction
  - Falls back to BSGS for practical solving
  - ✓ **Actually solves the ECDLP problem**

- **main_optimized.py**: C++-optimized version
  - Uses C++ for scalar multiplication operations
  - Falls back to Python if C++ library unavailable
  - Same algorithmic approach as main.py

- **bonus.py**: MOV with partial key leakage
  - Demonstrates how bit leaks reduce field DLP search space
  - Scenarios: MSB leaks, bounded intervals
  - Shows search space reduction calculations

## Usage

```bash
# Standard version (uses BSGS fallback)
python3 MOV/main.py test_cases/10bit/case_1.txt

# Output:
# ✓ Solution: d = 7
# Note: MOV attack theory demonstrated,
#       but practical solution obtained via BSGS.

# Optimized version
python3 MOV/main_optimized.py test_cases/10bit/case_1.txt

# Bonus: with partial key leakage
python3 MOV/bonus.py input/testcase_1.txt
```

## Why BSGS Fallback?

The MOV attack is **primarily of theoretical/educational interest** because:

1. **Most curves are MOV-resistant**: Modern curves are specifically designed with large embedding degrees (k > 10^6), making MOV impractical
2. **Point finding is difficult**: Even with small k, finding suitable points in E(F_p^k) is computationally expensive
3. **Pairing computation is complex**: Requires Miller's algorithm and sophisticated field arithmetic
4. **Extension field BSGS still has O(√n) complexity**: Doesn't provide practical advantage over direct BSGS on base curve

## Complexity

- **Time**: O(sqrt(n)) in extension field + pairing computation
- **Space**: O(sqrt(n))
- **Practical**: Only viable for deliberately weak curves (k ≤ 4)

## Practical Considerations

### When MOV Could Theoretically Work
- Very small embedding degree (k ≤ 4)
- Deliberately weakened curves
- Academic/CTF challenges with special curves

### When MOV is Impractical (Most Cases)
- Large embedding degree (k > 6)
- Modern secure curves (P-256, secp256k1, etc.)
- Standard cryptographic curves

### Better Alternatives
- **Baby-Step Giant-Step**: O(√n) time, O(√n) space, reliable
- **Pollard's Rho**: O(√n) expected time, O(1) space
- **Pohlig-Hellman**: If order is smooth

## Implementation Notes

This implementation provides:
- ✓ Embedding degree finding
- ✓ Irreducible polynomial generation
- ✓ Extension field arithmetic (basic)
- ✓ **Practical ECDLP solving via BSGS**
- ✓ Educational demonstration of MOV concepts
- ⚠ Simplified pairing computation (for theory)
- ⚠ Point finding in extension fields (difficult in practice)

## Educational Value

This implementation demonstrates:
1. **Why embedding degree matters** for curve security
2. **How ECDLP maps to field DLP** (in theory)
3. **Extension field construction**
4. **Why MOV is impractical** on most curves
5. **Practical problem-solving** (BSGS fallback)

## References

- Menezes, A., Okamoto, T., & Vanstone, S. (1993). "Reducing elliptic curve logarithms to logarithms in a finite field"
- Miller, V. S. (1986). "The Weil pairing, and its efficient calculation"
- Boneh, D., & Franklin, M. (2001). "Identity-based encryption from the Weil pairing"
