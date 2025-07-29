/**
 * CommonUtils.h - Shared utilities for number theory density computations
 * 
 * Contains common functionality used by both abundant and covering number 
 * density calculation programs, including prime utilities, mathematical 
 * constants, and shared algorithmic components.
 */

#ifndef COMMON_UTILS_H
#define COMMON_UTILS_H

#include <cfenv>
#include <cmath>
#include <string>
#pragma STDC FENV_ACCESS ON
#include <sstream>
#include <iostream>
#include <iomanip>
#include <gmpxx.h>
#include <mpfr.h>
#include <fstream>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <memory>
#include <vector>

// Algorithm configuration constants
#define RLEN 25                    // Length of lambda arrays for bounds computation
#define DEFAULT_PRINT_COUNT 20000000   // Default progress reporting interval

// Mathematical constants for number theory calculations
constexpr long double ABUNDANT_THRESHOLD = 2.0L;     // σ(n)/n ≥ 2 for abundant numbers (floating point)
constexpr long double COVERING_THRESHOLD = 2.0L;     // σ(n)/n ≥ 2 for covering numbers (floating point)
constexpr long double AK_BOUND_TARGET = 2.0L;        // Target value for AK bound computation
constexpr int ABUNDANT_THRESHOLD_INT = 2;            // σ(n)/n ≥ 2 for abundant numbers (exact arithmetic)

using std::cout;
using std::cerr;
using std::endl;
using std::string;

/**
 * Utility class for prime number operations
 * Provides efficient primality testing and prime generation
 */
class PrimeUtils {
public:
    /**
     * Check if a number is prime using trial division
     * Optimized for moderate-sized numbers with early exit conditions
     * @param num Number to test for primality
     * @return true if num is prime, false otherwise
     */
    static bool isPrime(long long num) {
        if (num <= 1) return false;
        if (num == 2) return true;
        if (num % 2 == 0) return false;
        
        // Trial division up to √num, checking only odd candidates
        for (long long i = 3; i * i <= num; i += 2) {
            if (num % i == 0) return false;
        }
        return true;
    }
    
    /**
     * Find the next prime greater than num
     * Uses incremental search with primality testing
     * @param num Starting number
     * @return Next prime after num
     */
    static long long nextPrime(long long num) {
        if (num < 2) return 2;
        
        // Start from next odd number if num is even
        num = (num % 2 == 0) ? num + 1 : num + 2;
        
        while (!isPrime(num)) {
            num += 2;  // Only check odd numbers
        }
        return num;
    }
};

/**
 * Mathematical utility functions for interval arithmetic and bounds
 */
class MathUtils {
public:
    /**
     * Set floating point rounding mode with error checking
     * @param mode Rounding mode (FE_UPWARD, FE_DOWNWARD, etc.)
     * @return true if successful, false otherwise
     */
    static bool setRoundingMode(int mode) {
        return std::fesetround(mode) == 0;
    }
    
    /**
     * Compute integer power efficiently using exponentiation by squaring
     * @param base Base value
     * @param exp Non-negative integer exponent
     * @return base^exp
     */
    template<typename T>
    static T intPower(T base, unsigned int exp) {
        T result = 1;
        while (exp > 0) {
            if (exp & 1) result *= base;
            base *= base;
            exp >>= 1;
        }
        return result;
    }
    
    /**
     * Safe conversion from mpz_class to built-in integer types
     * @param z GMP integer to convert
     * @return Converted value, or throws if out of range
     */
    static int safeMpzToInt(const mpz_class& z) {
        if (!z.fits_sint_p()) {
            throw std::overflow_error("mpz_class value too large for int conversion");
        }
        return static_cast<int>(z.get_si());
    }
};

/**
 * File I/O utilities for data persistence and recovery
 */
class FileUtils {
public:
    /**
     * Generate standardized filename for data files
     * @param prefix File prefix (e.g., "SavedData", "LamdasStore")
     * @param numPrimes Number of primes parameter
     * @param rlen R-length parameter
     * @param zpow Z-power parameter  
     * @param suffix File suffix (e.g., "out.dat", "in.dat")
     * @return Formatted filename string
     */
    static string generateDataFilename(const string& prefix, int numPrimes, 
                                      int rlen, int zpow, const string& suffix) {
        return prefix + std::to_string(numPrimes) + '-' + std::to_string(rlen) + 
               '-' + std::to_string(zpow) + suffix;
    }
    
    /**
     * Check if a file exists and is readable
     * @param filename Path to file
     * @return true if file exists and is accessible
     */
    static bool fileExists(const string& filename) {
        return std::filesystem::exists(filename) && 
               std::filesystem::is_regular_file(filename);
    }
    
    /**
     * Get file size in bytes
     * @param filename Path to file
     * @return File size, or 0 if file doesn't exist
     */
    static size_t getFileSize(const string& filename) {
        if (!fileExists(filename)) return 0;
        return std::filesystem::file_size(filename);
    }
};

/**
 * Progress reporting and logging utilities
 */
class ProgressUtils {
public:
    /**
     * Print progress message with timestamp
     * @param message Progress message to display
     * @param taskId Optional task identifier
     */
    static void logProgress(const string& message, int taskId = -1) {
        if (taskId >= 0) {
            cout << "[Task " << taskId << "] " << message << endl;
        } else {
            cout << message << endl;
        }
        cout.flush();
    }
    
    /**
     * Print computational statistics in standardized format
     * @param count Numbers processed
     * @param taskId Current task identifier
     * @param currentNum Current number being processed
     * @param primeIndex Current prime index
     */
    static void logComputationStats(unsigned long long count, int taskId, 
                                   const mpz_class& currentNum, int primeIndex) {
        cout << taskId << ' ' << currentNum << ' ' << count << ' ' << primeIndex << endl;
    }
    
    /**
     * Print final results in standardized format (legacy space-separated)
     * @param count Total numbers processed
     * @param lowerBound Lower bound estimate
     * @param upperBound Upper bound estimate  
     * @param w1Sum W₁ contribution
     * @param w2Sum W₂ contribution
     * @param densityAccounted Fraction of total density accounted for
     * @param givenUp Contribution from abandoned branches
     * @param lowerBoundUp Upper version of lower bound
     */
    static void printFinalResults(unsigned long long count, double lowerBound, 
                                 double upperBound, double w1Sum, double w2Sum,
                                 double densityAccounted, double givenUp, 
                                 double lowerBoundUp) {
        cout << std::setprecision(30) 
             << count << ' ' << lowerBound << ' ' << upperBound << ' ' 
             << w1Sum << ' ' << w2Sum << ' ' << densityAccounted << ' ' 
             << givenUp << ' ' << lowerBoundUp << endl;
    }
    
    /**
     * Print final results in human-readable format for publication
     * @param numberType Type of numbers being analyzed ("Abundant" or "Covering")
     * @param count Total numbers processed
     * @param lowerBound Lower bound estimate
     * @param upperBound Upper bound estimate  
     * @param w1Sum W₁ contribution
     * @param w2Sum W₂ contribution
     * @param densityAccounted Fraction of total density accounted for
     * @param givenUp Contribution from abandoned branches
     * @param lowerBoundUp Upper version of lower bound
     */
    static void printReadableResults(const string& numberType, unsigned long long count, 
                                   double lowerBound, double upperBound, double w1Sum, 
                                   double w2Sum, double densityAccounted, double givenUp, 
                                   double lowerBoundUp) {
        cout << std::fixed << std::setprecision(15);
        
        cout << "\n" << std::string(70, '=') << endl;
        cout << "FINAL RESULTS: " << numberType << " Number Density Bounds" << endl;
        cout << std::string(70, '=') << endl;
        
        cout << "\nCOMPUTATION STATISTICS:" << endl;
        cout << "  Numbers processed: " << std::setw(20) << count << endl;
        cout << "  Total density accounted: " << std::setw(15) << densityAccounted << endl;
        
        cout << "\nDENSITY BOUNDS:" << endl;
        cout << "  Lower bound:       " << std::setw(18) << lowerBound << endl;
        cout << "  Upper bound:       " << std::setw(18) << upperBound << endl;
        cout << "  Bound width:       " << std::setw(18) << (upperBound - lowerBound) << endl;
        cout << "  Relative precision: " << std::setw(17) << (upperBound - lowerBound) / lowerBound << endl;
        
        cout << "\nCONTRIBUTION BREAKDOWN:" << endl;
        cout << "  W₁ sum (early term): " << std::setw(15) << w1Sum << endl;
        cout << "  W₂ sum (power term): " << std::setw(15) << w2Sum << endl;
        cout << "  Given up branches:   " << std::setw(15) << givenUp << endl;
        cout << "  Lower bound (upper): " << std::setw(15) << lowerBoundUp << endl;
        
        cout << std::string(70, '=') << endl;
        
        // Also print legacy format for compatibility
        cout << "\nLegacy format (space-separated for scripts):" << endl;
        printFinalResults(count, lowerBound, upperBound, w1Sum, w2Sum, densityAccounted, givenUp, lowerBoundUp);
    }
};

/**
 * Memory management utilities for large array operations
 */
class MemoryUtils {
public:
    /**
     * Allocate and initialize 2D array of long doubles
     * @param rows Number of rows
     * @param cols Number of columns  
     * @return Unique pointer to allocated array
     */
    static std::unique_ptr<long double*[]> allocate2DArray(int rows, int cols) {
        auto array = std::make_unique<long double*[]>(rows);
        for (int i = 0; i < rows; i++) {
            array[i] = new long double[cols]();  // Zero-initialize
        }
        return array;
    }
    
    /**
     * Deallocate 2D array allocated by allocate2DArray
     * @param array Array to deallocate
     * @param rows Number of rows in array
     */
    static void deallocate2DArray(std::unique_ptr<long double*[]>& array, int rows) {
        if (array) {
            for (int i = 0; i < rows; i++) {
                delete[] array[i];
            }
            array.reset();
        }
    }
};

/**
 * Validation utilities for computational correctness
 */
class ValidationUtils {
public:
    /**
     * Validate that total density is approximately 1.0
     * @param totalDensity Computed total density
     * @param tolerance Acceptable tolerance (default 1e-10)
     * @return true if density is valid
     */
    static bool validateTotalDensity(long double totalDensity, 
                                    long double tolerance = 1e-10L) {
        return std::abs(totalDensity - 1.0L) <= tolerance;
    }
    
    /**
     * Validate computational parameters before starting
     * @param zpow Search depth parameter
     * @param numPrimes Number of primes
     * @param prec MPFR precision
     * @return true if parameters are valid
     */
    static bool validateParameters(int zpow, int numPrimes, int prec) {
        if (zpow < 10 || zpow > 100) {
            cerr << "Error: zpow must be between 10 and 100" << endl;
            return false;
        }
        if (numPrimes < 1 || numPrimes > 1000000) {
            cerr << "Error: numPrimes must be between 1 and 1000000" << endl;
            return false;
        }
        if (prec < 53 || prec > 1000) {
            cerr << "Error: precision must be between 53 and 1000 bits" << endl;
            return false;
        }
        return true;
    }
    
    /**
     * Check for integer overflow in computations
     * @param value Value to check
     * @param maxValue Maximum allowed value
     * @param context Description of where overflow occurred
     * @return true if value is within bounds
     */
    template<typename T>
    static bool checkOverflow(T value, T maxValue, const string& context) {
        if (value > maxValue) {
            cerr << "Error: Integer overflow in " << context 
                 << " (value=" << value << ", max=" << maxValue << ")" << endl;
            return false;
        }
        return true;
    }
};

#endif // COMMON_UTILS_H
