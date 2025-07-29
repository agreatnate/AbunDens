/**
 * DensityMathFramework.h - Shared mathematical framework for density calculations
 * 
 * This file contains the common mathematical infrastructure used by both abundant
 * and covering number density calculations, including:
 * - Prime sequence generation and density computation
 * - Lambda bound calculations using Deleglise's method
 * - AK bound computations for algorithmic termination
 * - High-precision interval arithmetic utilities
 * 
 * By centralizing these functions, we ensure consistency between the two programs
 * and make maintenance much easier.
 */

#ifndef DENSITY_MATH_FRAMEWORK_H
#define DENSITY_MATH_FRAMEWORK_H

#include "CommonUtils.h"

/**
 * Core mathematical data and operations for number theory density calculations
 * This class manages the shared computational infrastructure needed by both
 * abundant and covering number density programs.
 */
class DensityMathFramework {
public:
    // Configuration parameters
    int numPrimes;                        // Number of primes to consider
    int rbd;                             // R-bound parameter for lambda arrays
    long double z;                       // Search depth parameter (2^zpow)
    long double zz;                      // zpow³ for distributed computation
    int maxPrime;                        // Largest prime in consideration
    int prec;                           // MPFR precision in bits

    // Precomputed arrays for efficient computation
    std::unique_ptr<int[]> ourPrimes_;                    // Prime sequence
    std::unique_ptr<long double[]> ourDensitiesDup_;      // Density upper bounds
    std::unique_ptr<long double[]> ourDensitiesDdn_;      // Density lower bounds
    std::unique_ptr<int[]> piX_;                          // Prime counting function π(x)
    std::unique_ptr<long double*[]> lambdas_;             // Lambda bound arrays
    std::unique_ptr<long double*[]> lambdasRatio_;        // Lambda ratio arrays

    DensityMathFramework() = default;
    ~DensityMathFramework();

    // Non-copyable for safety (large arrays and unique ownership)
    DensityMathFramework(const DensityMathFramework&) = delete;
    DensityMathFramework& operator=(const DensityMathFramework&) = delete;

    /**
     * Initialize the complete mathematical framework
     * @param zpow Search depth parameter
     * @param numPrimes Number of consecutive primes to precompute
     * @param prec MPFR precision for high-accuracy computation
     * @param problemType String identifier for file naming ("Abund" or "Cov")
     */
    void initialize(int zpow, int numPrimes, int prec, const string& problemType);

    /**
     * Initialize prime list and density computation arrays
     * @param numPrimes Number of consecutive primes to precompute
     * @param prec MPFR precision for density computations
     */
    void initializePrimes(int numPrimes, int prec);
    
    /**
     * Initialize lambda bounds (load from file or compute fresh)
     * @param numPrimes Number of primes for bounds computation
     * @param prec MPFR precision for high-accuracy computation
     * @return Highest valid R-bound index
     */
    int initializeLambdas(int numPrimes, int prec);

    /**
     * Compute AK bounds using lambda values (critical performance path)
     * This is the core bound computation used throughout both algorithms
     * @param goal Target bound value (typically 2.0)
     * @param currentValueUP Current value upper bound (abundancy or covering index)
     * @param lambdas Lambda bound array for current prime
     * @param lambdasRatio Lambda ratio array for current prime
     * @return Computed bound value
     * @note Assumes FE_DOWNWARD rounding mode is already set by caller
     */
    inline long double akBound(long double goal, long double currentValueUP, 
                              const long double* lambdas, const long double* lambdasRatio) const;

    /**
     * Get density bounds for a given prime index
     * @param k Prime index
     * @return Pair of (lower_bound, upper_bound) for density
     */
    std::pair<long double, long double> getDensityBounds(int k) const {
        if (k < 0 || k > numPrimes) {
            throw std::out_of_range("Prime index out of range: " + std::to_string(k));
        }
        return {ourDensitiesDdn_[k], ourDensitiesDup_[k]};
    }

    /**
     * Get prime value by index
     * @param k Prime index (0 = special value 1, 1 = first prime 2, etc.)
     * @return Prime value
     */
    int getPrime(int k) const {
        if (k < 0 || k > numPrimes) {
            throw std::out_of_range("Prime index out of range: " + std::to_string(k));
        }
        return ourPrimes_[k];
    }

    /**
     * Get prime index for a given value using inverse pi function
     * @param value Integer value to look up
     * @return Prime index, or -1 if value exceeds maxPrime
     */
    int getPrimeIndex(int value) const {
        if (value < 0 || value > maxPrime) return -1;
        return piX_[value];
    }

    /**
     * Get lambda bounds for a given prime index
     * @param k Prime index
     * @return Pair of pointers to (lambdas, lambdasRatio) arrays
     */
    std::pair<const long double*, const long double*> getLambdaBounds(int k) const {
        if (k < 0 || k > numPrimes) {
            throw std::out_of_range("Prime index out of range: " + std::to_string(k));
        }
        return {lambdas_[k], lambdasRatio_[k]};
    }

    /**
     * Validate that the framework is properly initialized
     * @return true if all components are initialized
     */
    bool isInitialized() const {
        return ourPrimes_ && ourDensitiesDup_ && ourDensitiesDdn_ && 
               piX_ && lambdas_ && lambdasRatio_ && numPrimes > 0;
    }

private:
    /**
     * Compute Deleglise bounds for lambda values
     * @param curP Current prime for bound computation
     * @param Lambda Array of lambda values to update
     * @param prec MPFR precision for calculations
     */
    void delegliseBound(const mpz_class& curP, mpfr_t* Lambda, int prec) const;
    
    /**
     * Compute Euler multiplicand bounds for lambda computation
     * @param invP Reciprocal of current prime
     * @param denom Denominator for bound computation
     * @param Lambda Lambda array to update
     * @param curP Current prime
     * @param r Power parameter
     * @param rpow Power index
     * @param prec MPFR precision for calculations
     */
    void eulerMultiplicandBound(const mpq_class& invP, const mpz_class& denom, 
                               mpfr_t* Lambda, const mpz_class& curP, 
                               unsigned long long r, int rpow, int prec) const;
};

/**
 * Template class for shared computational variables and operations
 * This provides a common interface for accumulating results across both
 * abundant and covering number computations.
 * 
 * @tparam T The specific variable type (SmoothNumbersAbundVars or SmoothNumbersCovVars)
 */
template<typename T>
class DensityComputationBase {
protected:
    const DensityMathFramework& framework_;
    
public:
    explicit DensityComputationBase(const DensityMathFramework& framework) 
        : framework_(framework) {}

    /**
     * Handle cases where we've exhausted our prime set
     * @param invSmoothUp Inverse of smooth number (upper bound)
     * @param currentValueUP Current value upper bound (abundancy or covering index)
     * @param vars Accumulation variables
     */
    void reachPrimeBoundary(long double invSmoothUp, long double currentValueUP, T& vars) const {
        static const mpz_class dummySmoothNum(1);
        incrementCount(1, dummySmoothNum, framework_.numPrimes + 1, vars);

        vars.totalDensityAccountedUP += framework_.ourDensitiesDup_[framework_.numPrimes] * invSmoothUp;

        MathUtils::setRoundingMode(FE_DOWNWARD);
        long double akReturn = framework_.akBound(AK_BOUND_TARGET, currentValueUP, 
                                                 framework_.lambdas_[framework_.numPrimes], 
                                                 framework_.lambdasRatio_[framework_.numPrimes]);
        // Note: akBound sets rounding mode to FE_UPWARD before returning
        vars.givenUp += akReturn * framework_.ourDensitiesDup_[framework_.numPrimes] * invSmoothUp;
    }

    /**
     * Update counting statistics and progress reporting
     * @param inc Count increment
     * @param smoothNum Current smooth number being processed
     * @param k Current prime index
     * @param vars Accumulation variables
     */
    void incrementCount(int inc, const mpz_class& smoothNum, int k, T& vars) const {
        vars.count += inc;
        if (vars.count > vars.printCount) {
            ProgressUtils::logComputationStats(vars.count, vars.taskIndex, smoothNum, k);
            vars.printCount += 2000000;  // Update next reporting threshold
        }
    }

    /**
     * Add contribution to W₁ sum (early termination case)
     * @param smoothNum Current smooth number
     * @param invSmoothUp Inverse smooth number (upper bound)
     * @param akReturn AK bound return value
     * @param k Prime index
     * @param vars Accumulation variables
     */
    void addToW1Sum(const mpz_class& smoothNum, long double invSmoothUp, 
                   long double akReturn, int k, T& vars) const {
        incrementCount(1, smoothNum, k, vars);
        vars.totalDensityAccountedUP += framework_.ourDensitiesDup_[k] * invSmoothUp;
        vars.w1SumUP += akReturn * framework_.ourDensitiesDup_[k] * invSmoothUp;
    }

    /**
     * Check if computation should terminate at current prime level
     * Uses density bounds at k-1 and compares to k
     * @param currentValueUP Current value upper bound
     * @param k Prime index
     * @param invSmoothDown Inverse smooth number (lower bound)
     * @return true if computation should terminate early
     */
    bool shouldTerminateAtPrimeLevel(long double currentValueUP, int k, long double invSmoothDown) const {
        if (k <= 0 || k > framework_.numPrimes) return false;
        
        MathUtils::setRoundingMode(FE_DOWNWARD);
        long double akReturn = framework_.akBound(AK_BOUND_TARGET, currentValueUP, 
                                                 framework_.lambdas_[k - 1], 
                                                 framework_.lambdasRatio_[k - 1]);
        // Note: akBound sets rounding mode to FE_UPWARD before returning
        //MathUtils::setRoundingMode(FE_UPWARD);
        
        return (akReturn * framework_.ourDensitiesDdn_[k - 1] * framework_.z * invSmoothDown < k);
    }

    /**
     * Check if should give up on current prime power iteration
     * Uses density bounds at k and compares to k+1
     * @param currentValueUP Current value upper bound  
     * @param k Prime index
     * @param invSmoothDown Inverse smooth number (lower bound)
     * @param useDistributedBound Whether to use distributed computation bound (zz vs z)
     * @return true if should give up on this prime power
     */
    bool shouldGiveUpOnPrimePower(long double currentValueUP, int k, long double invSmoothDown, 
                                 bool useDistributedBound = false) const {
        if (k <= 0 || k > framework_.numPrimes) return false;
        
        MathUtils::setRoundingMode(FE_DOWNWARD);
        long double akReturn = framework_.akBound(AK_BOUND_TARGET, currentValueUP, 
                                                 framework_.lambdas_[k], 
                                                 framework_.lambdasRatio_[k]);
        // Note: akBound sets rounding mode to FE_UPWARD before returning
        
        long double bound = useDistributedBound ? framework_.zz : framework_.z;
        return (akReturn * framework_.ourDensitiesDdn_[k] * bound * invSmoothDown < k + 1);
    }
};

// ============================================================================
// Inline Implementation
// ============================================================================

inline long double DensityMathFramework::akBound(const long double goal, 
                                                 const long double currentValueUP, 
                                                 const long double* lambdas, 
                                                 const long double* lambdasRatio) const {
    /**
     * Critical performance path: compute AK bounds using lambda arrays
     * This function is called very frequently, so it's heavily optimized
     * 
     * Algorithm:
     * 1. Compute xd = goal / currentValueUP
     * 2. Check if we can't do any better than trivial (lambdas[0] > xd - 1)
     * 3. Iterate through lambda ratios, squaring xd at each step
     * 4. Return appropriate bound based on iteration results
     * 
     * NOTE: Assumes FE_DOWNWARD rounding mode is already set by caller **********
     */
    

    //Lambdas[r] stores values of  mu_{2^(r)} -1, the 2^(r)th moment minus 1
    //LambdasRatio[r] stores values of (mu_{2^(r+1)}-1) / (mu_{2^r} -1) - 1
    //So long as x^(2^r) is greater than this value we keep proceeding to the next value of r.

    long double xd = goal / currentValueUP;

    // Early termination check
    if (lambdas[0] > (xd - 1.0L)) {
        MathUtils::setRoundingMode(FE_UPWARD);
        return 1.0L;
    }

    // Main iteration loop
    int r;
    for (r = 1; r <= RLEN - 1; r++) {
        if (xd <= lambdasRatio[r]) {
            MathUtils::setRoundingMode(FE_UPWARD);
            return lambdas[r - 1] / (xd - 1.0L);
        }
        xd = xd * xd;  // Square for next iteration
    }
    
    // Final case: exhausted all ratios
    MathUtils::setRoundingMode(FE_UPWARD);
    return lambdas[r - 1] / (xd - 1.0L);
}

#endif // DENSITY_MATH_FRAMEWORK_H
