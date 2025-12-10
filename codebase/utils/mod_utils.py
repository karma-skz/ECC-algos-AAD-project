"""
Modular arithmetic utilities for ECC algorithms.
"""

from typing import List, Tuple


def extended_gcd(a: int, b: int) -> Tuple[int, int, int]:
    """
    Extended Euclidean Algorithm.
    
    Returns (gcd, x, y) such that a*x + b*y = gcd(a, b).
    """
    if b == 0:
        return (abs(a), 1 if a >= 0 else -1, 0)
    
    gcd, x1, y1 = extended_gcd(b, a % b)
    x = y1
    y = x1 - (a // b) * y1
    return (gcd, x, y)


def mod_inv(a: int, m: int) -> int:
    """
    Compute modular inverse of a modulo m.
    
    Args:
        a: The number to find inverse of
        m: The modulus (should be prime for guaranteed inverse)
    
    Returns:
        The modular inverse of a (mod m)
    
    Raises:
        ZeroDivisionError: If inverse doesn't exist
    """
    a %= m
    if a == 0:
        raise ZeroDivisionError(f"No inverse for 0 modulo {m}")
    
    gcd, x, _ = extended_gcd(a, m)
    if gcd != 1:
        raise ZeroDivisionError(f"No inverse: gcd({a}, {m}) = {gcd}")
    
    return x % m


def crt_combine(congruences: List[Tuple[int, int]]) -> Tuple[int, int]:
    """
    Chinese Remainder Theorem: combine congruences.
    
    Args:
        congruences: List of (remainder, modulus) pairs
    
    Returns:
        (solution, combined_modulus) where solution is the combined remainder
    
    Assumes moduli are pairwise coprime.
    """
    if not congruences:
        return (0, 1)
    
    result = 0
    product = 1
    
    for _, mod in congruences:
        product *= mod
    
    for remainder, mod in congruences:
        partial_product = product // mod
        _, inverse, _ = extended_gcd(partial_product, mod)
        result = (result + remainder * partial_product * inverse) % product
    
    return (result % product, product)
