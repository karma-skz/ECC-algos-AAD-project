"""
Bonus utilities for partial key leakage scenarios.
Simulates various types of information leaks and helper functions.
"""

import random
from typing import Tuple, List, Optional


class KeyLeakage:
    """Generate and manage partial key leakage scenarios."""
    
    @staticmethod
    def leak_lsb_bits(d: int, num_bits: int) -> Tuple[int, int]:
        """
        Leak the least significant bits of d.
        
        Args:
            d: Secret value
            num_bits: Number of LSB bits to leak
        
        Returns:
            (leaked_value, mask) where mask has num_bits set
        """
        mask = (1 << num_bits) - 1
        leaked = d & mask
        return leaked, mask
    
    @staticmethod
    def leak_msb_bits(d: int, num_bits: int, total_bits: int) -> Tuple[int, int]:
        """
        Leak the most significant bits of d.
        
        Args:
            d: Secret value
            num_bits: Number of MSB bits to leak
            total_bits: Total bit length of d
        
        Returns:
            (leaked_value, shift_amount)
        """
        # Handle edge case where secret is smaller than expected
        actual_bits = d.bit_length()
        if actual_bits < num_bits:
            # If secret is too small, just return the secret itself
            return d, 0
        
        shift = max(0, total_bits - num_bits)
        leaked = d >> shift
        return leaked, shift
    
    @staticmethod
    def leak_bounded_interval(d: int, n: int, leak_percentage: float = 0.1) -> Tuple[int, int]:
        """
        Leak that d is within a bounded interval.
        
        Args:
            d: Secret value
            n: Maximum value (order)
            leak_percentage: Percentage of n to use as interval size (0.0 to 1.0)
        
        Returns:
            (lower_bound, upper_bound)
        """
        interval_size = max(1, int(n * leak_percentage))
        half_interval = interval_size // 2
        
        lower = max(0, d - half_interval)
        upper = min(n - 1, d + half_interval)
        
        return lower, upper
    
    @staticmethod
    def leak_residues(d: int, moduli: List[int]) -> List[Tuple[int, int]]:
        """
        Leak residues of d modulo small primes.
        
        Args:
            d: Secret value
            moduli: List of moduli to leak residues for
        
        Returns:
            List of (residue, modulus) pairs
        """
        return [(d % m, m) for m in moduli]
    
    @staticmethod
    def leak_approximate_value(d: int, n: int, error_percentage: float = 0.05) -> Tuple[int, int]:
        """
        Leak an approximate value of d with some error.
        
        Args:
            d: Secret value
            n: Maximum value
            error_percentage: Maximum error as percentage of n
        
        Returns:
            (approximate_value, max_error)
        """
        max_error = max(1, int(n * error_percentage))
        error = random.randint(-max_error, max_error)
        approx = (d + error) % n
        return approx, max_error


def candidates_from_lsb_leak(leaked_lsb: int, mask: int, n: int) -> List[int]:
    """
    Generate candidate values consistent with LSB leak.
    
    Args:
        leaked_lsb: Known LSB bits
        mask: Bit mask for known bits
        n: Maximum value
    
    Returns:
        List of candidate values
    """
    candidates = []
    step = mask + 1  # Number of values between candidates
    
    for d in range(leaked_lsb, n, step):
        candidates.append(d)
    
    return candidates


def candidates_from_msb_leak(leaked_msb: int, shift: int, n: int) -> Tuple[int, int]:
    """
    Get range of candidates consistent with MSB leak.
    
    Args:
        leaked_msb: Known MSB bits
        shift: Number of bits shifted right
        n: Maximum value
    
    Returns:
        (lower_bound, upper_bound)
    """
    lower = leaked_msb << shift
    upper = min(n - 1, ((leaked_msb + 1) << shift) - 1)
    return lower, upper


def candidates_from_interval(lower: int, upper: int) -> Tuple[int, int]:
    """
    Get bounded interval for search.
    
    Args:
        lower: Lower bound
        upper: Upper bound
    
    Returns:
        (lower, upper) validated bounds
    """
    return max(0, lower), upper


def format_leak_info(leak_type: str, **kwargs) -> str:
    """
    Format leak information for display.
    
    Args:
        leak_type: Type of leak
        **kwargs: Leak-specific parameters
    
    Returns:
        Formatted string describing the leak
    """
    if leak_type == "lsb":
        num_bits = kwargs.get('num_bits', 0)
        leaked = kwargs.get('leaked', 0)
        return f"LSB Leak: {num_bits} bits known (value: {leaked:#b})"
    
    elif leak_type == "msb":
        num_bits = kwargs.get('num_bits', 0)
        leaked = kwargs.get('leaked', 0)
        return f"MSB Leak: {num_bits} bits known (value: {leaked:#b})"
    
    elif leak_type == "interval":
        lower = kwargs.get('lower', 0)
        upper = kwargs.get('upper', 0)
        size = upper - lower + 1
        return f"Bounded Interval: [{lower}, {upper}] (size: {size:,})"
    
    elif leak_type == "residues":
        residues = kwargs.get('residues', [])
        res_str = ", ".join(f"d≡{r} (mod {m})" for r, m in residues)
        return f"Residue Leak: {res_str}"
    
    elif leak_type == "approximate":
        approx = kwargs.get('approx', 0)
        error = kwargs.get('error', 0)
        return f"Approximate Value: {approx} ± {error}"
    
    else:
        return f"Unknown leak type: {leak_type}"


def calculate_speedup(original_time: float, leak_time: float) -> float:
    """Calculate speedup factor."""
    if leak_time == 0:
        return float('inf')
    return original_time / leak_time


def calculate_search_reduction(original_space: int, reduced_space: int) -> Tuple[float, str]:
    """
    Calculate search space reduction.
    
    Returns:
        (reduction_factor, percentage_string)
    """
    if reduced_space == 0:
        return float('inf'), "100.00%"
    
    reduction = original_space / reduced_space
    percentage = 100.0 * (1 - reduced_space / original_space)
    
    return reduction, f"{percentage:.2f}%"
