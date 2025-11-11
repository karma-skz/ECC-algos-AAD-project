#include <cstdint>

extern "C" {

// Extended GCD for modular inverse
int64_t mod_inv(int64_t a, int64_t m) {
    int64_t m0 = m, x0 = 0, x1 = 1;
    
    if (m == 1) return 0;
    
    while (a > 1) {
        int64_t q = a / m;
        int64_t t = m;
        m = a % m;
        a = t;
        t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }
    
    if (x1 < 0) x1 += m0;
    return x1;
}

// Modular multiplication to avoid overflow
int64_t mod_mult(int64_t a, int64_t b, int64_t p) {
    return ((__int128)a * b) % p;
}

// Point addition on elliptic curve
// Returns 0 if result is infinity, 1 if valid point
int point_add(int64_t x1, int64_t y1, int64_t x2, int64_t y2,
              int64_t a, int64_t p,
              int64_t* x3, int64_t* y3) {
    
    // Identity cases
    if (x1 < 0) {  // First point is infinity
        *x3 = x2;
        *y3 = y2;
        return (x2 >= 0) ? 1 : 0;
    }
    if (x2 < 0) {  // Second point is infinity
        *x3 = x1;
        *y3 = y1;
        return 1;
    }
    
    int64_t slope;
    
    if (x1 == x2) {
        if (y1 == y2) {
            // Point doubling
            int64_t numerator = mod_mult(3, mod_mult(x1, x1, p), p);
            numerator = (numerator + a) % p;
            int64_t denominator = mod_mult(2, y1, p);
            slope = mod_mult(numerator, mod_inv(denominator, p), p);
        } else {
            // Points are inverses
            return 0;  // Result is infinity
        }
    } else {
        // Point addition
        int64_t numerator = (y2 - y1 + p) % p;
        int64_t denominator = (x2 - x1 + p) % p;
        slope = mod_mult(numerator, mod_inv(denominator, p), p);
    }
    
    *x3 = (mod_mult(slope, slope, p) - x1 - x2 + 2*p) % p;
    *y3 = (mod_mult(slope, (x1 - *x3 + p), p) - y1 + p) % p;
    
    return 1;
}

// Scalar multiplication using double-and-add
// Returns 0 if result is infinity, 1 if valid point
int scalar_mult(int64_t k, int64_t gx, int64_t gy,
                int64_t a, int64_t b, int64_t p,
                int64_t* result_x, int64_t* result_y) {
    
    if (k == 0) {
        return 0;  // Infinity
    }
    
    // Initialize result to infinity
    int64_t rx = -1, ry = -1;
    int result_is_infinity = 1;
    
    // Current point (used for doubling)
    int64_t px = gx, py = gy;
    
    while (k > 0) {
        if (k & 1) {
            // Add current point to result
            if (result_is_infinity) {
                rx = px;
                ry = py;
                result_is_infinity = 0;
            } else {
                int64_t new_x, new_y;
                int valid = point_add(rx, ry, px, py, a, p, &new_x, &new_y);
                if (!valid) {
                    return 0;  // Result is infinity
                }
                rx = new_x;
                ry = new_y;
            }
        }
        
        // Double the current point
        if (k > 1) {
            int64_t new_x, new_y;
            int valid = point_add(px, py, px, py, a, p, &new_x, &new_y);
            if (!valid) {
                // Doubling resulted in infinity (shouldn't happen for valid curve)
                return 0;
            }
            px = new_x;
            py = new_y;
        }
        
        k >>= 1;
    }
    
    if (result_is_infinity) {
        return 0;
    }
    
    *result_x = rx;
    *result_y = ry;
    return 1;
}

}  // extern "C"
