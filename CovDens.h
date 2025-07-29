/**
 * CovDens19.h - Covering number density calculation using shared framework
 * 
 * This simplified version leverages the DensityMathFramework for all common
 * mathematical operations, focusing only on covering-number-specific logic.
 * Code size reduced by ~50% while maintaining identical functionality.
 */

#ifndef COVDENS19_H
#define COVDENS19_H

#include "DensityMathFramework.h"

/**
 * Structure representing a divisor with its multiplicity information
 * Used for Bell number computations in covering index calculations
 */
struct Divisor {
    mpz_class d;      // The divisor value
    int omega;        // Number of distinct prime factors of d
    
    Divisor(const mpz_class& di, int o) : d(di), omega(o) {}
};

/**
 * Container for variables used during covering number traversal
 * Mirrors the abundant version but with covering-specific naming
 */
class SmoothNumbersCovVars {
public:
    // Lower bounds for covering number density (interval arithmetic)
    long double lowBdCovering;      // Lower bound with downward rounding
    long double lowBdCoveringUp;    // Lower bound with upward rounding
    
    // Weighted sums for different algorithmic cases
    long double w1SumUP;            // W₁ sum (early termination cases)
    long double w2SumUP;            // W₂ sum (deep search cases)
    
    // Tracking and validation variables
    long double totalDensityAccountedUP;  // Total density processed (should → 1)
    long double givenUp;                  // Contribution from abandoned branches
    unsigned long long count;             // Numbers processed
    unsigned long long printCount;        // Progress printing threshold
    unsigned long long givenUpcount;      // Count of abandoned tasks
    int taskIndex;                        // Current parallel task ID

    // Reusable GMP objects to minimize allocation overhead
    mpq_class tempMaxCoveringQ;
    mpz_class tempMaxCoveringZ;

    SmoothNumbersCovVars();
    SmoothNumbersCovVars(const SmoothNumbersCovVars& other);
    SmoothNumbersCovVars& operator=(const SmoothNumbersCovVars& other);
    
    void mergeWith(const SmoothNumbersCovVars& other);
};

/**
 * Enhanced mathematical framework for covering numbers
 * Extends the base framework with Bell number computations
 */
class CoveringMathFramework : public DensityMathFramework {
private:
    // Bell numbers cache for covering index computation
    std::vector<std::vector<long long>> bellCache_;

public:
    /**
     * Initialize covering-specific components including Bell numbers
     * @param zpow Search depth parameter
     * @param numPrimes Number of consecutive primes to precompute
     * @param prec MPFR precision for high-accuracy computation
     */
    void initializeCovering(int zpow, int numPrimes, int prec);

    /**
     * Retrieve Bell number B(j,k) from cache with bounds checking
     * @param j First parameter of Bell number
     * @param k Second parameter of Bell number
     * @return Bell number value, or throws exception if out of range
     */
    long long getBellNumber(int j, int k) const;

private:
    void initializeBellCache();
};

/**
 * Covering number traversal using shared mathematical framework
 * This class inherits common operations and focuses on covering-specific logic
 */
class TraverseCoveringNumbers : public DensityComputationBase<SmoothNumbersCovVars> {
private:
    CoveringMathFramework framework_;
    
public:
    TraverseCoveringNumbers() : DensityComputationBase(framework_) {}
    
    /**
     * Main entry point: compute density bounds for covering numbers
     * @param zpow Search depth parameter (actual depth = 2^zpow)
     * @param numPrimes Number of consecutive primes to consider
     * @param prec MPFR precision in bits for high-accuracy computation
     */
    void calculate(int zpow, int numPrimes, int prec);
    
private:
    // === Initialization and Control Methods ===
    int determineTaskCount(SmoothNumbersCovVars& final);
    void loadExistingResults(const string& infilename, const string& outfilename,
                            std::vector<bool>& processedTasks, SmoothNumbersCovVars& final,
                            std::ofstream& fout);
    void computeRemainingTasks(int numTasks, const std::vector<bool>& processedTasks,
                              SmoothNumbersCovVars& final, std::ofstream& fout);
    void outputFinalResults(const SmoothNumbersCovVars& final) const;
    
    // === Covering-Specific Mathematical Functions ===
    
    /**
     * Process a newly discovered covering number
     * @param invSmoothUp Inverse smooth number (upper bound)
     * @param invSmoothDown Inverse smooth number (lower bound)
     * @param k Prime index
     * @param vars Accumulation variables
     */
    inline void newCoveringNumber(long double invSmoothUp, long double invSmoothDown, 
                          int k, SmoothNumbersCovVars& vars) const;
    
    /**
     * Process an interval of covering numbers
     * @param invSmoothUp Inverse smooth number (upper bound)
     * @param invSmoothDown Inverse smooth number (lower bound)
     * @param initk Initial prime index
     * @param k Final prime index
     * @param vars Accumulation variables
     */
    inline void newCoveringInterval(long double invSmoothUp, long double invSmoothDown, 
                            int initk, int k, SmoothNumbersCovVars& vars) const;

    /**
     * Compute covering index using Bell number formula
     * Formula: 1 + (ℓ-1)/ℓ + (1/ℓ) * Σ_{d|b, d>1} B(ω(d), τ(ℓ))/d
     * @param smoothNum Current smooth number
     * @param sunDivCount Number of divisors of current smooth part
     * @param ell Parameter ℓ in covering index formula
     * @param divisors List of divisors with multiplicity information
     * @return Computed covering index
     */
    long double computeCoveringIndexDirect(const mpz_class& smoothNum, int sunDivCount, 
                                          const mpz_class& ell, const std::vector<Divisor>& divisors) const;
    
    /**
     * Compute shifted covering index sum for optimization
     * @param sunDivCount Number of divisors of current smooth part
     * @param divisors List of divisors with multiplicity information
     * @return Shifted covering index contribution
     */
    inline long double computeShiftedCoveringIndexSum(int sunDivCount, 
                                              const std::vector<Divisor>& divisors) const;

    /**
     * Handle automatic covering detection for small primes
     * @param k Current prime index (will be updated)
     * @param sunDivCount Number of divisors
     * @param invSmoothUp Inverse smooth number (upper bound)
     * @param invSmoothDown Inverse smooth number (lower bound)
     * @param vars Accumulation variables
     * @param isDistrib Whether this is distributed computation
     * @return Updated prime index k
     */
    int handleAutomaticCovering(int k, int sunDivCount, long double invSmoothUp, 
                               long double invSmoothDown, SmoothNumbersCovVars& vars, 
                               bool isDistrib) const;

    // === Core Traversal Methods ===
    
    /**
     * Main recursive traversal for covering number enumeration
     * Uses shared framework for all common mathematical operations
     * @param smoothNum Current smooth number being built
     * @param invSmoothUp Inverse of smooth number (upper bound)
     * @param invSmoothDown Inverse of smooth number (lower bound)
     * @param k Current prime index
     * @param sunDivCount Number of divisors in smooth part
     * @param coveringIndex Current covering index value
     * @param ell Parameter ℓ for covering computation
     * @param divisors Current divisor list with multiplicities
     * @param vars Result accumulation variables
     */
    void smoothNumbersCov(const mpz_class& smoothNum, long double invSmoothUp, 
                         long double invSmoothDown, int k, int sunDivCount, 
                         long double coveringIndex, const mpz_class& ell, 
                         const std::vector<Divisor>& divisors, SmoothNumbersCovVars& vars);
    
    /**
     * Distributed version for parallel task decomposition
     * Uses shared framework for all common mathematical operations
     * @param smoothNum Current smooth number being built
     * @param invSmoothUp Inverse of smooth number (upper bound)
     * @param invSmoothDown Inverse of smooth number (lower bound)
     * @param k Current prime index
     * @param sunDivCount Number of divisors in smooth part
     * @param ell Parameter ℓ for covering computation
     * @param divisors Current divisor list with multiplicities
     * @param vars Result accumulation variables
     */
    void smoothNumbersCovDistrib(const mpz_class& smoothNum, long double invSmoothUp, 
                                long double invSmoothDown, int k, int sunDivCount, 
                                const mpz_class& ell, const std::vector<Divisor>& divisors, 
                                SmoothNumbersCovVars& vars);
};

#endif // COVDENS19_H
