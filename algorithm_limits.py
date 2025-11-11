#!/usr/bin/env python3
"""
Algorithm Feasibility Configuration

Defines which algorithms can handle which bit sizes of ECDLP problems.
"""

# Algorithm feasibility limits (in bits)
ALGORITHM_LIMITS = {
    'BruteForce': {
        'practical': 14,  # Practically feasible
        'maximum': 40,    # Theoretically can test, but very slow
        'description': 'Exhaustive search. Practical up to 14 bits, theoretical up to 40 bits.'
    },
    'BabyStep': {
        'practical': 40,  # Memory-dependent but generally feasible
        'maximum': 40,
        'description': 'Baby-step Giant-step. Feasible up to 40 bits (memory-dependent).'
    },
    'PohligHellman': {
        'practical': 40,  # Depends on factorization
        'maximum': 40,
        'description': 'Pohlig-Hellman. Feasible up to 40 bits (factorization-dependent).'
    },
    'PollardRho': {
        'practical': 16,  # Probabilistic, works sometimes
        'maximum': 20,    # Beyond this: NOT FEASIBLE
        'description': 'Pollard Rho. Probabilistic, NOT FEASIBLE beyond 20 bits.'
    },
    'LasVegas': {
        'practical': 14,  # Probabilistic, works sometimes
        'maximum': 20,    # Beyond this: NOT FEASIBLE
        'description': 'Las Vegas. Probabilistic, NOT FEASIBLE beyond 20 bits.'
    }
}

def is_feasible(algorithm, bits):
    """Check if an algorithm is feasible for a given bit size."""
    if algorithm not in ALGORITHM_LIMITS:
        return False
    return bits <= ALGORITHM_LIMITS[algorithm]['maximum']

def is_practical(algorithm, bits):
    """Check if an algorithm is practically feasible for a given bit size."""
    if algorithm not in ALGORITHM_LIMITS:
        return False
    return bits <= ALGORITHM_LIMITS[algorithm]['practical']

def get_feasible_range(algorithm):
    """Get the feasible bit range for an algorithm."""
    if algorithm not in ALGORITHM_LIMITS:
        return None
    return (10, ALGORITHM_LIMITS[algorithm]['maximum'])

def get_practical_range(algorithm):
    """Get the practical bit range for an algorithm."""
    if algorithm not in ALGORITHM_LIMITS:
        return None
    return (10, ALGORITHM_LIMITS[algorithm]['practical'])

def print_algorithm_limits():
    """Print a summary of algorithm limits."""
    print("=" * 80)
    print("ECDLP Algorithm Feasibility Limits")
    print("=" * 80)
    
    for algo, limits in ALGORITHM_LIMITS.items():
        print(f"\n{algo}:")
        print(f"  Practical: Up to {limits['practical']} bits")
        print(f"  Maximum:   Up to {limits['maximum']} bits")
        print(f"  {limits['description']}")
    
    print("\n" + "=" * 80)
    print("Test Cases: 10-40 bits (31 sizes × 5 cases = 155 total)")
    print("=" * 80)

if __name__ == "__main__":
    print_algorithm_limits()
