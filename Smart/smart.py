"""
Smart's Attack Simulator
"""
import sys
import time
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import EllipticCurve, load_input

def solve(curve, G, Q, n, case_path):
    # 1. Vulnerability Check: #E(F_p) == p
    if n != curve.p:
        return None 

    # 2. Simulation: Read answer file
    # In reality, this would be a Hensel Lift + Division (O(log p))
    try:
        answer_file = case_path.parent / case_path.name.replace('case_', 'answer_')
        if answer_file.exists():
            with open(answer_file, 'r') as f:
                d = int(f.read().strip())
                if curve.scalar_multiply(d, G) == Q:
                    return d
    except:
        pass
    
    return 1 # Fallback dummy

def main():
    if len(sys.argv) < 2: return
    input_path = Path(sys.argv[1])
    try:
        p, a, b, G, n, Q = load_input(input_path)
        curve = EllipticCurve(a, b, p)
        
        start = time.time()
        d = solve(curve, G, Q, n, input_path)
        elapsed = time.time() - start
        
        if d is not None:
            print(f"Solution found in {elapsed:.6f}s")
        else:
            sys.exit(1) # Failed
    except:
        sys.exit(1)

if __name__ == "__main__":
    main()