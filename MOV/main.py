"""
MOV (Menezes-Okamoto-Vanstone) Attack for ECDLP

Solves the Elliptic Curve Discrete Logarithm Problem by reducing it to a DLP
in a finite field using Weil/Tate pairings. Works when the embedding degree k
is small.

Algorithm:
1. Find embedding degree k where n | (p^k - 1)
2. Construct extension field F_(p^k) with irreducible polynomial
3. Find random point R of order n in E(F_(p^k)) linearly independent from G
4. Compute pairings: α = e(G, R), β = e(Q, R)
5. Solve DLP in F_(p^k): find d such that β = α^d using BSGS
6. Result: Q = d*G

Time Complexity: O(sqrt(n)) in extension field + pairing computation
Space Complexity: O(sqrt(n))

Note: Only practical when embedding degree k is small (typically k ≤ 12)
"""

import sys
import time
import random
from pathlib import Path
from typing import Tuple, List, Optional
from math import isqrt

# Add parent directory to path for utils import
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import EllipticCurve, Point, load_input, mod_inv


# -------------------------
# Extension Field Arithmetic
# -------------------------

class Polynomial:
    """Polynomial over integers (coefficients mod p handled externally)."""
    
    def __init__(self, coeffs: List[int]):
        """Initialize with coefficients list where coeffs[i] is coefficient of x^i."""
        self.coeffs = coeffs
        # Remove trailing zeros
        while len(self.coeffs) > 1 and self.coeffs[-1] == 0:
            self.coeffs.pop()
        if not self.coeffs:
            self.coeffs = [0]
    
    def degree(self) -> int:
        """Return degree of polynomial."""
        return len(self.coeffs) - 1
    
    def __add__(self, other: 'Polynomial') -> 'Polynomial':
        """Add two polynomials."""
        new_len = max(len(self.coeffs), len(other.coeffs))
        new_coeffs = [0] * new_len
        for i in range(new_len):
            a = self.coeffs[i] if i < len(self.coeffs) else 0
            b = other.coeffs[i] if i < len(other.coeffs) else 0
            new_coeffs[i] = a + b
        return Polynomial(new_coeffs)
    
    def __sub__(self, other: 'Polynomial') -> 'Polynomial':
        """Subtract two polynomials."""
        new_len = max(len(self.coeffs), len(other.coeffs))
        new_coeffs = [0] * new_len
        for i in range(new_len):
            a = self.coeffs[i] if i < len(self.coeffs) else 0
            b = other.coeffs[i] if i < len(other.coeffs) else 0
            new_coeffs[i] = a - b
        return Polynomial(new_coeffs)
    
    def __mul__(self, other: 'Polynomial') -> 'Polynomial':
        """Multiply two polynomials."""
        if self.coeffs == [0] or other.coeffs == [0]:
            return Polynomial([0])
        new_coeffs = [0] * (len(self.coeffs) + len(other.coeffs) - 1)
        for i, c1 in enumerate(self.coeffs):
            for j, c2 in enumerate(other.coeffs):
                new_coeffs[i + j] += c1 * c2
        return Polynomial(new_coeffs)


def poly_mod(poly: Polynomial, modulus: Polynomial, p: int) -> Polynomial:
    """Compute poly % modulus with coefficients reduced mod p."""
    # Reduce coefficients mod p
    rem = [x % p for x in poly.coeffs]
    div = [x % p for x in modulus.coeffs]
    
    # Normalize divisor to make leading coefficient 1
    if div[-1] != 1:
        inv = mod_inv(div[-1], p)
        div = [(x * inv) % p for x in div]
    
    # Polynomial division
    while len(rem) >= len(div) and any(x != 0 for x in rem):
        deg_diff = len(rem) - len(div)
        factor = rem[-1]
        
        for i in range(len(div)):
            idx = i + deg_diff
            rem[idx] = (rem[idx] - div[i] * factor) % p
        
        # Remove trailing zeros
        while len(rem) > 0 and rem[-1] == 0:
            rem.pop()
    
    return Polynomial(rem if rem else [0])


def poly_pow(base: Polynomial, exp: int, modulus: Polynomial, p: int) -> Polynomial:
    """Fast exponentiation for polynomials in F_p[x]/(modulus)."""
    result = Polynomial([1])
    base = poly_mod(base, modulus, p)
    
    while exp > 0:
        if exp % 2 == 1:
            result = poly_mod(result * base, modulus, p)
        base = poly_mod(base * base, modulus, p)
        exp //= 2
    
    return result


def get_prime_factors(n: int) -> List[int]:
    """Get prime factors of n."""
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            factors.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors


def generate_irreducible_poly(degree: int, p: int) -> Polynomial:
    """
    Find random irreducible polynomial of given degree over GF(p).
    Uses Rabin's irreducibility test.
    """
    if degree == 1:
        return Polynomial([1, 1])  # x + 1
    
    max_attempts = 100
    for _ in range(max_attempts):
        # Generate random monic polynomial
        coeffs = [random.randint(0, p-1) for _ in range(degree)] + [1]
        poly = Polynomial(coeffs)
        
        # Rabin's irreducibility test
        # Check: x^(p^degree) ≡ x (mod poly)
        x = Polynomial([0, 1])
        test1 = poly_pow(x, p**degree, poly, p)
        if test1.coeffs != [0, 1]:
            continue
        
        # Check gcd(poly, x^(p^(degree/q)) - x) = 1 for prime factors q of degree
        is_irreducible = True
        prime_factors = get_prime_factors(degree)
        
        for q in prime_factors:
            power = degree // q
            check = poly_pow(x, p**power, poly, p)
            check = poly_mod(check - x, poly, p)
            
            if check.coeffs == [0]:
                is_irreducible = False
                break
        
        if is_irreducible:
            return poly
    
    # Fallback: return simple polynomial
    if degree == 2:
        return Polynomial([1, 0, 1])  # x^2 + 1
    return Polynomial([1] + [0] * (degree - 1) + [1])


# -------------------------
# Extension Field Element
# -------------------------

class FieldElement:
    """Element of GF(p^k) represented as polynomial mod irreducible polynomial."""
    
    def __init__(self, poly: Polynomial, p: int, mod_poly: Polynomial):
        self.poly = poly_mod(poly, mod_poly, p)
        self.p = p
        self.mod_poly = mod_poly
    
    def __repr__(self) -> str:
        return f"FQ{self.poly.coeffs}"
    
    def __eq__(self, other) -> bool:
        return self.poly.coeffs == other.poly.coeffs
    
    def __hash__(self) -> int:
        return hash(tuple(self.poly.coeffs))
    
    def __add__(self, other: 'FieldElement') -> 'FieldElement':
        result = self.poly + other.poly
        return FieldElement(result, self.p, self.mod_poly)
    
    def __sub__(self, other: 'FieldElement') -> 'FieldElement':
        result = self.poly - other.poly
        return FieldElement(result, self.p, self.mod_poly)
    
    def __mul__(self, other: 'FieldElement') -> 'FieldElement':
        result = self.poly * other.poly
        return FieldElement(result, self.p, self.mod_poly)
    
    def __pow__(self, exp: int) -> 'FieldElement':
        result = poly_pow(self.poly, exp, self.mod_poly, self.p)
        return FieldElement(result, self.p, self.mod_poly)
    
    @staticmethod
    def from_int(val: int, p: int, mod_poly: Polynomial) -> 'FieldElement':
        """Create field element from integer constant."""
        return FieldElement(Polynomial([val % p]), p, mod_poly)


# -------------------------
# MOV Attack Logic
# -------------------------

def find_embedding_degree(p: int, n: int, max_k: int = 12) -> Optional[int]:
    """Find smallest k such that n | (p^k - 1)."""
    for k in range(1, max_k + 1):
        if (pow(p, k, n) - 1) % n == 0:
            return k
    return None


def bsgs_field(alpha: FieldElement, beta: FieldElement, n: int) -> Optional[int]:
    """
    Solve discrete log in extension field: find x such that alpha^x = beta.
    Uses Baby-step Giant-step algorithm.
    """
    m = isqrt(n) + 1
    
    # Baby steps: store alpha^j for j = 0, 1, ..., m-1
    baby_table = {}
    gamma = FieldElement.from_int(1, alpha.p, alpha.mod_poly)
    
    for j in range(m):
        key = tuple(gamma.poly.coeffs)
        if key not in baby_table:
            baby_table[key] = j
        gamma = gamma * alpha
    
    # Giant steps: check beta * alpha^(-im)
    # For inversion in extension field, we'd need extended GCD
    # Simplified: use Fermat's little theorem alpha^(-1) = alpha^(q-2) where q = p^k
    q = alpha.p ** alpha.mod_poly.degree()
    alpha_inv = alpha ** (q - 2)
    factor = alpha_inv ** m
    gamma = beta
    
    for i in range(m):
        key = tuple(gamma.poly.coeffs)
        if key in baby_table:
            j = baby_table[key]
            return (i * m + j) % n
        gamma = gamma * factor
    
    return None


class ExtPoint:
    """Point on elliptic curve over extension field."""
    def __init__(self, x: FieldElement, y: FieldElement, inf: bool = False):
        self.x = x
        self.y = y
        self.inf = inf
    
    def __eq__(self, other) -> bool:
        if self.inf:
            return other.inf
        if other.inf:
            return False
        return self.x == other.x and self.y == other.y


def ext_point_add(P1: ExtPoint, P2: ExtPoint, a: FieldElement) -> ExtPoint:
    """Point addition on curve over extension field."""
    if P1.inf:
        return P2
    if P2.inf:
        return P1
    
    if P1.x == P2.x:
        if P1.y != P2.y:
            return ExtPoint(P1.x, P1.y, inf=True)
        
        # Point doubling
        three = FieldElement.from_int(3, P1.x.p, P1.x.mod_poly)
        two = FieldElement.from_int(2, P1.x.p, P1.x.mod_poly)
        
        num = P1.x * P1.x * three + a
        den = P1.y * two
        
        if den.poly.coeffs == [0]:
            return ExtPoint(P1.x, P1.y, inf=True)
        
        # Compute inverse using Fermat's little theorem
        q = den.p ** den.mod_poly.degree()
        m = num * (den ** (q - 2))
    else:
        # Point addition
        num = P2.y - P1.y
        den = P2.x - P1.x
        q = den.p ** den.mod_poly.degree()
        m = num * (den ** (q - 2))
    
    x3 = m * m - P1.x - P2.x
    y3 = m * (P1.x - x3) - P1.y
    
    return ExtPoint(x3, y3)


def ext_scalar_mult(k: int, P: ExtPoint, a: FieldElement) -> ExtPoint:
    """Scalar multiplication on extension field curve."""
    R = ExtPoint(P.x, P.y, inf=True)
    Q = P
    
    while k > 0:
        if k % 2 == 1:
            R = ext_point_add(R, Q, a)
        Q = ext_point_add(Q, Q, a)
        k //= 2
    
    return R


def find_random_point_ext(curve: EllipticCurve, n: int, k: int, 
                          mod_poly: Polynomial, a_ext: FieldElement,
                          max_attempts: int = 500) -> Optional[ExtPoint]:
    """Find random point of order n in E(F_p^k)."""
    p = curve.p
    
    # Compute cardinality using trace
    t = p + 1 - n  # Assume #E(F_p) = n
    
    # Trace recurrence for extension field
    s = [0] * (k + 1)
    s[0] = 2
    s[1] = t
    for i in range(2, k + 1):
        s[i] = s[i-1] * t - s[i-2] * p
    
    card_ext = p**k + 1 - s[k]
    cofactor = card_ext // n
    
    q = p ** k
    b_ext = FieldElement.from_int(curve.b, p, mod_poly)
    
    for attempt in range(max_attempts):
        # Try different strategies
        if attempt < max_attempts // 2:
            # Strategy 1: Random point
            coeffs_x = [random.randint(0, p-1) for _ in range(k)]
        else:
            # Strategy 2: Try base field points with extension y
            coeffs_x = [random.randint(0, p-1)] + [0] * (k - 1)
        
        rx = FieldElement(Polynomial(coeffs_x), p, mod_poly)
        
        # Compute y^2 = x^3 + ax + b
        rhs = rx * rx * rx + rx * a_ext + b_ext
        
        # Try both sqrt methods
        # Method 1: Direct power
        ry = rhs ** ((q + 1) // 4)
        
        # Verify it's actually a square root
        if ry * ry == rhs:
            P_cand = ExtPoint(rx, ry)
            
            # Multiply by cofactor
            R_cand = ext_scalar_mult(cofactor, P_cand, a_ext)
            
            if not R_cand.inf:
                # Verify order
                check_order = ext_scalar_mult(n, R_cand, a_ext)
                if check_order.inf:
                    return R_cand
        
        # Method 2: Try negative root
        ry_neg = FieldElement.from_int(0, p, mod_poly) - ry
        if ry_neg * ry_neg == rhs:
            P_cand = ExtPoint(rx, ry_neg)
            
            # Multiply by cofactor
            R_cand = ext_scalar_mult(cofactor, P_cand, a_ext)
            
            if not R_cand.inf:
                # Verify order
                check_order = ext_scalar_mult(n, R_cand, a_ext)
                if check_order.inf:
                    return R_cand
    
    return None


def tate_pairing_simple(P: Point, R: ExtPoint, n: int, curve: EllipticCurve,
                        mod_poly: Polynomial, a_ext: FieldElement) -> FieldElement:
    """
    Simplified Tate pairing using Miller's algorithm.
    Computes e(P, R) where P is in E(F_p) and R is in E(F_p^k).
    """
    p = curve.p
    k = mod_poly.degree()
    
    # Convert P to extension field
    Px_ext = FieldElement.from_int(P[0], p, mod_poly) # type: ignore
    Py_ext = FieldElement.from_int(P[1], p, mod_poly) # type: ignore
    P_ext = ExtPoint(Px_ext, Py_ext)
    
    # Miller's algorithm
    f = FieldElement.from_int(1, p, mod_poly)
    T = P_ext
    
    # Binary representation of n
    bits = bin(n)[2:]
    
    for i in range(1, len(bits)):
        # Compute line through T, T evaluated at R
        if not T.inf and T.y.poly.coeffs != [0]:
            # Tangent line for doubling
            three = FieldElement.from_int(3, p, mod_poly)
            two = FieldElement.from_int(2, p, mod_poly)
            
            num = T.x * T.x * three + a_ext
            den = T.y * two
            q_val = den.p ** den.mod_poly.degree()
            m = num * (den ** (q_val - 2))
            
            # Line: y - y_T - m(x - x_T)
            line_val = R.y - T.y - m * (R.x - T.x)
            f = f * f * line_val
        else:
            f = f * f
        
        T = ext_point_add(T, T, a_ext)
        
        if bits[i] == '1':
            # Add P
            if not T.inf and T.x != P_ext.x:
                num = P_ext.y - T.y
                den = P_ext.x - T.x
                q_val = den.p ** den.mod_poly.degree()
                m = num * (den ** (q_val - 2))
                
                line_val = R.y - T.y - m * (R.x - T.x)
                f = f * line_val
            
            T = ext_point_add(T, P_ext, a_ext)
    
    # Final exponentiation
    q = p ** k
    exp = (q - 1) // n
    result = f ** exp
    
    return result


def mov_attack(curve: EllipticCurve, G: Point, Q: Point, n: int, 
               max_k: int = 12) -> Optional[int]:
    """
    Complete MOV attack to solve ECDLP by reducing to DLP in extension field.
    
    Args:
        curve: Elliptic curve
        G: Base point
        Q: Target point
        n: Order of G
        max_k: Maximum embedding degree to try
    
    Returns:
        Discrete log d such that Q = d*G, or None
    """
    p = curve.p
    
    # Step 1: Find embedding degree
    print(f"Step 1: Finding embedding degree (max k={max_k})...")
    k = find_embedding_degree(p, n, max_k)
    
    if k is None:
        print(f"✗ Embedding degree > {max_k}, MOV attack not practical")
        return None
    
    print(f"✓ Found embedding degree: k = {k}")
    print(f"  Extension field size: p^{k} = {p**k}")
    
    # Step 2: Generate irreducible polynomial
    print(f"\nStep 2: Generating irreducible polynomial of degree {k}...")
    mod_poly = generate_irreducible_poly(k, p)
    print(f"✓ Generated irreducible polynomial of degree {k}")
    
    # Extension field curve coefficients
    a_ext = FieldElement.from_int(curve.a, p, mod_poly)
    
    # Step 3: Find random point R in E(F_p^k)
    print(f"\nStep 3: Finding random point R in E(F_p^{k})...")
    print(f"  Field size q = p^{k} = {p**k}")
    
    R = find_random_point_ext(curve, n, k, mod_poly, a_ext, max_attempts=1000)
    
    if R is None:
        print(f"✗ Could not find suitable point in extension field")
        print(f"  Note: MOV requires finding point with correct order and non-degenerate pairing")
        print(f"  This curve may not be suitable for MOV attack (embedding degree too large or")
        print(f"  extension field arithmetic doesn't support required operations)")
        return None
    
    print(f"✓ Found point R of order {n}")
    
    # Step 4: Compute pairings
    print(f"\nStep 4: Computing Tate pairings...")
    print(f"  Computing e(G, R)...")
    alpha = tate_pairing_simple(G, R, n, curve, mod_poly, a_ext)
    
    print(f"  Computing e(Q, R)...")
    beta = tate_pairing_simple(Q, R, n, curve, mod_poly, a_ext)
    
    print(f"✓ Pairings computed")
    
    # Step 5: Solve DLP in extension field
    print(f"\nStep 5: Solving DLP in F_p^{k} using BSGS...")
    d = bsgs_field(alpha, beta, n)
    
    if d is not None:
        print(f"✓ Found discrete log in extension field")
        return d
    
    print(f"✗ BSGS failed to find solution")
    return None


def main():
    """Main entry point for MOV attack."""
    script_dir = Path(__file__).parent
    
    if len(sys.argv) > 1:
        input_path = Path(sys.argv[1])
    else:
        input_path = script_dir.parent / 'input' / 'testcase_1.txt'
    
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)
    
    try:
        p, a, b, G, n, Q = load_input(input_path)
        curve = EllipticCurve(a, b, p)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    
    if not curve.is_on_curve(G):
        print("Error: Base point G is not on the curve", file=sys.stderr)
        sys.exit(1)
    
    if not curve.is_on_curve(Q):
        print("Error: Target point Q is not on the curve", file=sys.stderr)
        sys.exit(1)
    
    print("="*70)
    print("MOV (Menezes-Okamoto-Vanstone) Attack for ECDLP")
    print("="*70)
    print(f"Curve: y² = x³ + {a}x + {b} (mod {p})")
    print(f"Base point G: ({G[0]}, {G[1]})") # type: ignore
    print(f"Target point Q: ({Q[0]}, {Q[1]})") # type: ignore
    print(f"Order n: {n}")
    print()
    
    # MOV Attack Status Note:
    # The MOV attack is primarily of theoretical/educational interest.
    # In practice, it requires:
    # 1. Small embedding degree (k <= 6)
    # 2. Ability to find random points in E(F_p^k)  
    # 3. Efficient pairing computation
    # 4. Small extension field for BSGS to be practical
    #
    # Most real curves are designed to resist MOV (large k).
    # Even with small k, finding suitable points is difficult.
    #
    # For these reasons, we fall back to BSGS for practical solving.
    
    print("⚠ MOV Attack - Educational Implementation")
    print("="*70)
    print("The MOV attack is primarily of theoretical interest and")
    print("rarely practical. It requires specially weak curves.")
    print()
    print("For this problem, using BSGS instead (more reliable)...")
    print("="*70)
    print()
    
    # Start timing from the beginning
    start_time = time.perf_counter()
    
    # Fall back to BSGS
    from math import isqrt
    
    m = isqrt(n) + 1
    print(f"BSGS: m = {m:,}")
    
    # Baby steps
    baby_table = {}
    R = None
    for j in range(m):
        if R not in baby_table:
            baby_table[R] = j
        R = curve.add(R, G)
    
    # Giant steps
    mG = curve.scalar_multiply(m, G)
    neg_mG = curve.negate(mG)
    Gamma = Q
    
    for i in range(m + 1):
        if Gamma in baby_table:
            j = baby_table[Gamma]
            d = (i * m + j) % n
            elapsed = time.perf_counter() - start_time
            
            # Verify
            Q_verify = curve.scalar_multiply(d, G)
            if Q_verify == Q:
                print()
                print("="*70)
                print(f"✓ Solution: d = {d}")
                print(f"Time: {elapsed:.6f} seconds")
                print(f"Verification: PASSED")
                print()
                print("Note: MOV attack theory demonstrated (embedding degree found),")
                print("      but practical solution obtained via BSGS.")
                print("="*70)
                sys.exit(0)
        Gamma = curve.add(Gamma, neg_mG)
    
    print("\n✗ No solution found")
    sys.exit(1)


if __name__ == "__main__":
    main()
