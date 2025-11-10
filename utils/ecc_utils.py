"""
Elliptic Curve utilities and operations.
"""

from typing import Optional, Tuple
from .mod_utils import mod_inv

Point = Optional[Tuple[int, int]]  # None represents the point at infinity


class EllipticCurve:
    """
    Represents an elliptic curve y^2 = x^3 + ax + b (mod p) over a prime field F_p.
    
    Attributes:
        a: Curve coefficient a
        b: Curve coefficient b
        p: Prime modulus defining the field
    """
    
    def __init__(self, a: int, b: int, p: int):
        """
        Initialize an elliptic curve.
        
        Args:
            a: Coefficient a in the curve equation
            b: Coefficient b in the curve equation
            p: Prime modulus
        
        Raises:
            ValueError: If curve is singular (discriminant is zero)
        """
        self.a = a % p
        self.b = b % p
        self.p = p
        
        if not self.is_valid():
            raise ValueError(f"Singular curve: 4a^3 + 27b^2 ≡ 0 (mod {p})")
    
    def is_valid(self) -> bool:
        """Check if the curve is non-singular (discriminant ≠ 0)."""
        discriminant = (4 * pow(self.a, 3, self.p) + 27 * pow(self.b, 2, self.p)) % self.p
        return discriminant != 0
    
    def is_on_curve(self, P: Point) -> bool:
        """Check if point P is on the curve."""
        if P is None:  # Point at infinity
            return True
        
        x, y = P
        if not (0 <= x < self.p and 0 <= y < self.p):
            return False
        
        lhs = (y * y) % self.p
        rhs = (x * x * x + self.a * x + self.b) % self.p
        return lhs == rhs
    
    def negate(self, P: Point) -> Point:
        """Return the additive inverse of point P."""
        if P is None:
            return None
        x, y = P
        return (x % self.p, (-y) % self.p)
    
    def add(self, P: Point, Q: Point) -> Point:
        """
        Add two points P and Q on the curve using the group law.
        
        Args:
            P: First point
            Q: Second point
        
        Returns:
            P + Q on the curve
        """
        # Handle identity element
        if P is None:
            return Q
        if Q is None:
            return P
        
        x1, y1 = P
        x2, y2 = Q
        
        # Check if P + Q = O (point at infinity)
        if x1 == x2 and (y1 + y2) % self.p == 0:
            return None
        
        # Calculate slope
        if P != Q:
            # Point addition: λ = (y2 - y1) / (x2 - x1)
            numerator = (y2 - y1) % self.p
            denominator = (x2 - x1) % self.p
            slope = (numerator * mod_inv(denominator, self.p)) % self.p
        else:
            # Point doubling: λ = (3x1^2 + a) / (2y1)
            if y1 == 0:
                return None  # Tangent is vertical
            numerator = (3 * x1 * x1 + self.a) % self.p
            denominator = (2 * y1) % self.p
            slope = (numerator * mod_inv(denominator, self.p)) % self.p
        
        # Calculate result
        x3 = (slope * slope - x1 - x2) % self.p
        y3 = (slope * (x1 - x3) - y1) % self.p
        
        return (x3, y3)
    
    def scalar_multiply(self, k: int, P: Point) -> Point:
        """
        Compute k * P using double-and-add algorithm.
        
        Args:
            k: Scalar multiplier
            P: Point to multiply
        
        Returns:
            k * P on the curve
        """
        if P is None:
            return None
        
        if k == 0:
            return None
        
        if k < 0:
            # Handle negative scalars: (-k)P = -(kP)
            return self.negate(self.scalar_multiply(-k, P))
        
        # Double-and-add algorithm
        result = None
        addend = P
        
        while k > 0:
            if k & 1:
                result = self.add(result, addend)
            addend = self.add(addend, addend)
            k >>= 1
        
        return result
    
    def find_order(self, G: Point, max_iter: int = 100000) -> Optional[int]:
        """
        Find the order of point G by repeated addition.
        
        Args:
            G: Generator point
            max_iter: Maximum iterations to prevent infinite loops
        
        Returns:
            Order of G, or None if not found within max_iter
        """
        if G is None:
            return 1
        
        R = G
        for k in range(1, max_iter + 1):
            if R is None:
                return k
            R = self.add(R, G)
        
        return None
    
    def __repr__(self) -> str:
        return f"EllipticCurve(y^2 = x^3 + {self.a}x + {self.b} mod {self.p})"
