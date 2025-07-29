# Bounding the density of abundant numbers and covering numbers.

This repository contains high-performance C++ programs for computing rigorous density bounds of two specific classes of integers.

## Programs

### Abundant Number Density Calculator (`AbunDens`)
Computes rigorous upper and lower bounds for the density of abundant numbers. A positive integer n is **abundant** if σ(n)/n ≥ 2, where σ(n) is the sum of divisors function.

### Covering Number Density Calculator (`CovDens`)
Computes rigorous upper bounds for the density of covering numbers. A positive integer n is **covering** if a distinct covering system of the integers can be formed using the divisors of n greater than one as moduli.

## Features

- **Rigorous Interval Arithmetic**: Uses high-precision MPFR arithmetic with configurable precision for calculation of moments
- **Rigorous Interval Arithmetic**: Uses long doubles with consistent rounding modes set for density calculations 
- **Parallel Processing**: OpenMP-based parallelization for efficient computation on multi-core systems
- **Persistent Computation**: Automatic save/resume functionality for long-running calculations

## Prerequisites

- **C++17 compatible compiler** (GCC 7+ recommended)
- **GMP** (GNU Multiple Precision Arithmetic Library)
- **MPFR** (Multiple Precision Floating-Point Reliably)
- **OpenMP** for parallel computation

### Installing Dependencies

**Ubuntu/Debian:**
```bash
sudo apt-get install build-essential libgmp-dev libmpfr-dev
```

**macOS (via Homebrew):**
```bash
brew install gmp mpfr
```

**CentOS/RHEL:**
```bash
sudo yum install gcc-c++ gmp-devel mpfr-devel
```

## Building

The project uses a comprehensive Makefile with multiple build targets:

```bash
# Build both programs (optimized)
make

# Build individual programs
make abundant
make covering

## Usage

Both programs accept similar command-line arguments:

```bash
# Basic usage with default precision (80 bits)
./bin/abundens <zpow> <numPrimes>

# Custom precision
./bin/abundens <zpow> <numPrimes> [precision]

# Examples
./bin/abundens 15 100000        # Default 80-bit precision
./bin/abundens 15 100000 120    # 120-bit precision
```

### Parameters

- **zpow**: Search depth parameter 
  - Controls search depth as the parameter z=2^zpow
  - Higher values = more precision but slower computation
  
- **numPrimes**: Number of consecutive primes to consider 
  - More primes = better bounds but longer computation time
  
- **precision**: MPFR precision in bits (optional, default: 80, range: 53-1000)
  - 80 bits seems sufficent for current calculations.

## Mathematical Background

### Abundant Numbers
An integer n > 1 is abundant if σ(n)/n ≥ 2, where σ(n) = Σ(d|n) d is the sum of divisors. The density of abundant numbers among positive integers is approximately 0.2476...

### Covering Numbers
An integer n is covering if a distinct covering system of the integers can be formed using the divisors of n greater than one as moduli. A covering system is a collection of arithmetic progressions whose union covers all integers. The computational test uses covering index c'(n) ≥ 2, computed using generalize Bell numbers.

## Output Format

```
FINAL RESULTS: Abundant Number Density Bounds
=============================================

COMPUTATION STATISTICS:
  Numbers processed:            12345678
  Total density accounted:      1.00000000002

DENSITY BOUNDS:
  Lower bound:       0.2474xxxxxxxxxxxxxx
  Upper bound:       0.2477xxxxxxxxxxxxxx
  Bound width:       0.0003xxxxxxxxxxxxxx

CONTRIBUTION BREAKDOWN:
  W₁ sum (early term):   0.123456789012345
  W₂ sum (power term):   0.098765432109876
  Given up branches:     0.025206340112346
  Lower bound (upper):   0.247428561234568
```

Note the Total density accounted should exceed 1 (it is designed to round up at every calculation.)

## File Structure

```
├── AbunDens.cpp/.h          # Abundant number calculator
├── CovDens.cpp/.h           # Covering number calculator  
├── DensityMathFramework.cpp/.h  # Shared mathematical framework
├── CommonUtils.h            # Utility functions and constants
├── BellNums.h              # Precomputed Bell numbers for covering calculations
├── Makefile                # Build system
└── README.md               # This file
```


## Citation
Please see the related paper for details on these calculations.

If you use these programs in research, please cite:
```
Nathan McNew and Jai Setty. On the density of covering numbers and abundant numbers.. 
```

## License

This software is released under the GNU General Public License v3.0 (GPLv3). See the [LICENSE](LICENSE) file for full details.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.


