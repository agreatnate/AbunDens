/**
 * AbunDens19.h - Abundant number density calculation using shared framework
 * 
 * This simplified version leverages the DensityMathFramework for all common
 * mathematical operations, focusing only on abundant-number-specific logic.
 * Code size reduced by ~60% while maintaining identical functionality.
 */

#ifndef ABUNDENS19_H
#define ABUNDENS19_H

#include "DensityMathFramework.h"

/**
 * Container for variables used during abundant number traversal
 * Identical to previous version but inherits shared functionality
 */
class SmoothNumbersAbundVars {
public:
    // Lower bounds for abundant number density (interval arithmetic)
    long double lowBdAbundant;      // Lower bound with downward rounding
    long double lowBdAbundantUp;    // Lower bound with upward rounding
    
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
    mpq_class tempMaxAbundQ;
    mpz_class tempMaxAbundZ;

    SmoothNumbersAbundVars();
    SmoothNumbersAbundVars(const SmoothNumbersAbundVars& other);
    SmoothNumbersAbundVars& operator=(const SmoothNumbersAbundVars& other);
    
    void mergeWith(const SmoothNumbersAbundVars& other);
};

/**
 * Abundant number traversal using shared mathematical framework
 * This class inherits common operations and focuses on abundant-specific logic
 */
class TraverseAbundantNumbers : public DensityComputationBase<SmoothNumbersAbundVars> {
private:
    DensityMathFramework framework_;
    
public:
    TraverseAbundantNumbers() : DensityComputationBase(framework_) {}
    
    /**
     * Main entry point: compute density bounds for abundant numbers
     * @param zpow Search depth parameter (actual depth = 2^zpow)
     * @param numPrimes Number of consecutive primes to consider
     * @param prec MPFR precision in bits for high-accuracy computation
     */
    void calculate(int zpow, int numPrimes, int prec);
    
private:
    // === Initialization and Control Methods ===
    int determineTaskCount(SmoothNumbersAbundVars& final);
    void loadExistingResults(const string& infilename, const string& outfilename,
                            std::vector<bool>& processedTasks, SmoothNumbersAbundVars& final,
                            std::ofstream& fout);
    void computeRemainingTasks(int numTasks, const std::vector<bool>& processedTasks,
                              SmoothNumbersAbundVars& final, std::ofstream& fout);
    void outputFinalResults(const SmoothNumbersAbundVars& final) const;
    
    // === Abundant-Specific Mathematical Functions ===
    
    /**
     * Process a newly discovered abundant number
     * @param invSmoothUp Inverse smooth number (upper bound)
     * @param invSmoothDown Inverse smooth number (lower bound) 
     * @param k Prime index
     * @param vars Accumulation variables
     */
    inline void newAbundantNumber(long double invSmoothUp, long double invSmoothDown, 
                          int k, SmoothNumbersAbundVars& vars) const;
    
    /**
     * Process an interval of abundant numbers
     * @param invSmoothUp Inverse smooth number (upper bound)
     * @param invSmoothDown Inverse smooth number (lower bound)
     * @param initk Initial prime index
     * @param k Final prime index  
     * @param vars Accumulation variables
     */
    inline void newAbundantInterval(long double invSmoothUp, long double invSmoothDown, 
                            int initk, int k, SmoothNumbersAbundVars& vars) const;
    

    /**
     * Handle initial abundant number detection and interval processing
     * @param smoothNum Current smooth number
     * @param curAbundancyQ Current abundancy (exact rational)
     * @param k Current prime index (will be updated)
     * @param invSmoothUp Inverse smooth number (upper bound)
     * @param invSmoothDown Inverse smooth number (lower bound)
     * @param vars Accumulation variables
     * @param isDistrib Whether this is distributed computation
     * @return Updated prime index k
     */
    int handleAbundantDetection(const mpz_class& smoothNum, const mpq_class& curAbundancyQ, int k,
                               long double invSmoothUp, long double invSmoothDown, 
                               SmoothNumbersAbundVars& vars, bool isDistrib) const;

    // === Core Traversal Methods ===
    
    /**
     * Main recursive traversal for abundant number enumeration
     * Uses shared framework for all common mathematical operations
     * @param smoothNum Current smooth number being built
     * @param invSmoothUp Inverse of smooth number (upper bound)
     * @param invSmoothDown Inverse of smooth number (lower bound)
     * @param k Current prime index
     * @param curAbundancyUP Current abundancy upper bound
     * @param curAbundancyQ Current abundancy (exact rational)
     * @param vars Result accumulation variables
     */
    void smoothNumbersAbund(const mpz_class& smoothNum, long double invSmoothUp, 
                           long double invSmoothDown, int k, long double curAbundancyUP, 
                           const mpq_class& curAbundancyQ, SmoothNumbersAbundVars& vars);
    
    /**
     * Distributed version for parallel task decomposition
     * Uses shared framework for all common mathematical operations
     * @param smoothNum Current smooth number being built
     * @param invSmoothUp Inverse of smooth number (upper bound)
     * @param invSmoothDown Inverse of smooth number (lower bound)
     * @param k Current prime index
     * @param curAbundancyUP Current abundancy upper bound
     * @param curAbundancyQ Current abundancy (exact rational)
     * @param vars Result accumulation variables
     */
    void smoothNumbersAbundDistrib(const mpz_class& smoothNum, long double invSmoothUp, 
                                  long double invSmoothDown, int k, long double curAbundancyUP, 
                                  const mpq_class& curAbundancyQ, SmoothNumbersAbundVars& vars);
};

#endif // ABUNDENS19_H
