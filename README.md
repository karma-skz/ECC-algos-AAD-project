# Elliptic Curve Discrete Logarithm Problem (ECDLP) Solvers

A comprehensive, modular implementation of various algorithms for solving the Elliptic Curve Discrete Logarithm Problem. This project demonstrates different algorithmic approaches with clean, well-documented, and maintainable code.

## 🚀 Quick Start

```bash
# Test an algorithm
cd BruteForce
python3 main.py input/test_1.txt

# Compare all algorithms
python3 compare_algorithms.py

# Generate test cases
cd utils
python3 generator.py
```

## 📚 Algorithms Implemented

1. **Brute Force** - O(n) exhaustive search
2. **Baby-Step Giant-Step** - O(√n) time/space tradeoff  
3. **Pohlig-Hellman** - Exploits smooth order factorization
4. **Pollard Rho** - O(√n) probabilistic with O(1) space
5. **Las Vegas** - Linear algebra based probabilistic approach

## 🏗️ Project Structure

```
├── utils/              # Shared utilities (ECC operations, modular math, I/O)
├── BruteForce/        # Exhaustive search implementation
├── BabyStep/          # Baby-Step Giant-Step (BSGS)
├── PohligHellman/     # Pohlig-Hellman algorithm
├── PollardRho/        # Pollard's Rho with Floyd cycle detection
├── LasVegas/          # Probabilistic linear algebra approach
└── compare_algorithms.py  # Benchmark all algorithms
```

## ✨ Key Features

- **Modular Design**: Shared utilities eliminate code duplication
- **Clean Code**: Comprehensive docstrings, type hints, clear naming
- **Verified Results**: All solutions verified against Q = d·G
- **Professional Quality**: Production-ready code structure
- **Educational**: Clear algorithm descriptions and complexity analysis
- **No Dependencies**: Uses only Python standard library

## 📖 Documentation

See [PROJECT_README.md](PROJECT_README.md) for comprehensive documentation including:
- Detailed algorithm explanations
- API documentation
- Input file format
- Usage examples
- Performance comparisons

## 🧪 Testing

All algorithms have been tested and verified:

```bash
# Individual tests
cd BruteForce && python3 main.py input/test_1.txt
cd BabyStep && python3 main.py input/test_1.txt  
cd PohligHellman && python3 main.py input/test_1.txt
cd PollardRho && python3 main.py input/test_1.txt
cd LasVegas && python3 main.py input/test_1.txt

# All algorithms on same input
python3 compare_algorithms.py
```

## 📊 Sample Output

```
Solving ECDLP using Baby-Step Giant-Step...
Curve: y^2 = x^3 + 497x + 1768 (mod 9739)
G = (1804, 5368), Q = (3138, 1774), n = 9735
m = ceil(sqrt(n)) = 99

==================================================
Solution: d = 1234
Time: 0.000228 seconds
Verification: PASSED
==================================================
```

## 🎓 Academic Context

This project was developed for the Advanced Algorithm Design (AAD) course, demonstrating:
- Algorithm design and analysis
- Data structure optimization
- Cryptographic mathematics
- Software engineering best practices

## 📝 Input Format

All algorithms use a standardized 5-line format:

```
<p>              # Prime modulus
<a> <b>          # Curve coefficients
<Gx> <Gy>        # Base point
<n>              # Order
<Qx> <Qy>        # Target point
```

## 🛠️ Requirements

- Python 3.7+
- No external dependencies

## 👨‍💻 Code Quality

- Type hints throughout
- Comprehensive docstrings
- Proper error handling
- Consistent code style
- Efficient algorithms
- Memory-conscious design

## 📈 Performance

| Algorithm | Time | Space | Best For |
|-----------|------|-------|----------|
| Brute Force | O(n) | O(1) | n < 10⁴ |
| BSGS | O(√n) | O(√n) | n < 10¹² |
| Pohlig-Hellman | O(∑√qᵢ) | O(√max qᵢ) | Smooth n |
| Pollard Rho | O(√n) | O(1) | General |
| Las Vegas | Polynomial | O(n'²) | Research |

## 🔒 Security Note

These implementations are for **educational purposes only**. For production cryptographic applications, use established libraries like `cryptography` or `PyCryptodome`.

## 📄 License

Academic project - AAD Course, Semester 2-1

---

**Made with ❤️ for learning and understanding cryptographic algorithms**



### input format-
```
p
a b
Gx Gy
n
Qx Qy
```

### output is the secret key `d` such that Q = d*G

## TODOs:
- Implement Pohlig-Hellman
- Find a dataset and adapt the I/O likewise.
- Add theory of each algo
- Bonus for each algo

- Actually look at the code and see what value the "no of attempts" kind of variables should take to balance compute taken and the accuracy