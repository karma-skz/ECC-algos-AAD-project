#!/usr/bin/env python3
"""
Auto-Breaker: Vulnerability Scanner & Exploit Router
Automatically detects the weakness in an Elliptic Curve and launches the appropriate attack.
"""

import sys
import time
import subprocess
from pathlib import Path
from utils.io_utils import load_input

# --- Configuration ---
# Threshold for "Smoothness" (for Pohlig-Hellman)
# If the largest prime factor is <= this, we consider it PH-vulnerable.
SMOOTH_BOUND = 100000 

def print_banner():
    print(r"""
    _   _   _ _____ ___  ____  ____  _____    _    _  _______ ____  
   / \ | | | |_   _/ _ \| __ )|  _ \| ____|  / \  | |/ / ____|  _ \ 
  / _ \| | | | | || | | |  _ \| |_) |  _|   / _ \ | ' /|  _| | |_) |
 / ___ \ |_| | | || |_| | |_) |  _ <| |___ / ___ \| . \| |___|  _ < 
/_/   \_\___/  |_| \___/|____/|_| \_\_____/_/   \_\_|\_\_____|_| \_\
                                                                    
    [+] ECC Vulnerability Scanner & Exploit Router
    """)

def is_smooth(n, bound):
    """
    Checks if n is B-smooth (all prime factors <= bound).
    Returns True if smooth, False otherwise.
    """
    temp_n = n
    d = 2
    while d * d <= temp_n and d <= bound:
        while temp_n % d == 0:
            temp_n //= d
        d += 1
    
    # If temp_n > 1, the remaining factor is a prime.
    # If that prime is > bound, then it's not smooth.
    if temp_n > bound:
        return False
    return True

def analyze_and_attack(file_path):
    print(f"[*] Loading target: {file_path}")
    try:
        p, a, b, G, n, Q = load_input(Path(file_path))
    except Exception as e:
        print(f"[!] Error loading file: {e}")
        return

    print(f"[*] Curve Parameters:")
    print(f"    p = {p} ({p.bit_length()} bits)")
    print(f"    n = {n}")
    print(f"    a = {a}, b = {b}")
    print("-" * 50)
    print("[*] Scanning for vulnerabilities...")
    time.sleep(0.5) # Dramatic pause

    attack_script = None
    attack_name = None
    
    # 1. Check for Anomalous Curve (Smart's Attack)
    if n == p:
        print("\033[91m[!] VULNERABILITY DETECTED: Anomalous Curve (#E = p)\033[0m")
        print("    -> Attack: Smart's Attack (Linear Time)")
        attack_script = "Smart/main.py"
        attack_name = "Smart's Attack"

    # 2. Check for Supersingular Curve (MOV Attack)
    # Common case: y^2 = x^3 + b over F_p where p = 2 mod 3 (j-invariant 0)
    # Or simply check if trace is 0 (n = p + 1)
    is_supersingular_form = (p % 3 == 2) and (a == 0)
    if n == p + 1 or is_supersingular_form:
        print("\033[91m[!] VULNERABILITY DETECTED: Supersingular Curve (Trace = 0)\033[0m")
        print("    -> Attack: MOV Attack (Reduction to DLP in F_p^k)")
        attack_script = "MOV/main_optimized.py"
        attack_name = "MOV Attack"

    # 3. Check for Smooth Order (Pohlig-Hellman)
    elif is_smooth(n, SMOOTH_BOUND):
        print("\033[91m[!] VULNERABILITY DETECTED: Smooth Order\033[0m")
        print(f"    -> Attack: Pohlig-Hellman (Decomposition)")
        attack_script = "PohligHellman/main_optimized.py"
        attack_name = "Pohlig-Hellman"

    # 4. Fallback: Generic (Pollard's Rho)
    else:
        print("\033[93m[-] No specific structural weakness found.\033[0m")
        print("    -> Attack: Pollard's Rho (Generic / O(sqrt(n)))")
        attack_script = "PollardRho/main_optimized.py"
        attack_name = "Pollard's Rho"

    print("-" * 50)
    print(f"[*] Launching {attack_name}...")
    print("-" * 50)

    # Execute the attack
    start_time = time.perf_counter()
    try:
        # Use the same python interpreter
        cmd = [sys.executable, attack_script, str(file_path)]
        result = subprocess.run(cmd, capture_output=False, text=True)
        
        end_time = time.perf_counter()
        if result.returncode == 0:
            print("-" * 50)
            print(f"\033[92m[+] Attack Successful!\033[0m")
            print(f"    Total Time: {end_time - start_time:.4f}s")
        else:
            print(f"\033[91m[!] Attack Failed with return code {result.returncode}\033[0m")
            
    except Exception as e:
        print(f"[!] Execution Error: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python auto_break.py <test_case_file>")
        sys.exit(1)
    
    print_banner()
    analyze_and_attack(sys.argv[1])
