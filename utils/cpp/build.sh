#!/bin/bash
# Build script for common ECC C++ library
g++ -O3 -shared -fPIC ecc_fast.cpp -o ecc_fast.so
