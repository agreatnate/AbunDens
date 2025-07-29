/**
 * CovDens.cpp - Covering number density calculation
 * 
 * This file implements algorithms to compute rigorous upper bounds for the 
 * density of covering numbers. An integer n is covering if a distinct covering
 * system of the integers can be formed using the divisors of n greater than one
 * as moduli. The implementation uses smooth number enumeration, Bell number
 * computations, interval arithmetic, and parallel processing to achieve
 * high-precision density bounds for research in computational number theory.
 */

#include "CovDens.h"
#include "BellNums.h"
#include "DensityMathFramework.h"
#include "CommonUtils.h"
#include <omp.h>

// ============================================================================
// SmoothNumbersCovVars Implementation
// ============================================================================

SmoothNumbersCovVars::SmoothNumbersCovVars() 
    : lowBdCovering(0.0L), lowBdCoveringUp(0.0L), w1SumUP(0.0L), w2SumUP(0.0L),
      totalDensityAccountedUP(0.0L), givenUp(0.0L), count(0), 
      printCount(DEFAULT_PRINT_COUNT), givenUpcount(0), taskIndex(0) {
}

SmoothNumbersCovVars::SmoothNumbersCovVars(const SmoothNumbersCovVars& other) 
    : lowBdCovering(other.lowBdCovering), lowBdCoveringUp(other.lowBdCoveringUp),
      w1SumUP(other.w1SumUP), w2SumUP(other.w2SumUP),
      totalDensityAccountedUP(other.totalDensityAccountedUP),
      givenUp(other.givenUp), count(other.count), printCount(other.printCount),
      givenUpcount(other.givenUpcount), taskIndex(other.taskIndex),
      tempMaxCoveringQ(other.tempMaxCoveringQ), tempMaxCoveringZ(other.tempMaxCoveringZ) {
}

SmoothNumbersCovVars& SmoothNumbersCovVars::operator=(const SmoothNumbersCovVars& other) {
    if (this != &other) {
        lowBdCovering = other.lowBdCovering;
        lowBdCoveringUp = other.lowBdCoveringUp;
        w1SumUP = other.w1SumUP;
        w2SumUP = other.w2SumUP;
        totalDensityAccountedUP = other.totalDensityAccountedUP;
        givenUp = other.givenUp;
        count = other.count;
        printCount = other.printCount;
        givenUpcount = other.givenUpcount;
        taskIndex = other.taskIndex;
        tempMaxCoveringQ = other.tempMaxCoveringQ;
        tempMaxCoveringZ = other.tempMaxCoveringZ;
    }
    return *this;
}

void SmoothNumbersCovVars::mergeWith(const SmoothNumbersCovVars& other) {
    // Use rigorous interval arithmetic for bound accumulation
    // ROUNDING MODE VALIDATION: Lower bounds must use FE_DOWNWARD
    MathUtils::setRoundingMode(FE_DOWNWARD);
    lowBdCovering += other.lowBdCovering;
    
    // ROUNDING MODE VALIDATION: Upper bounds must use FE_UPWARD
    MathUtils::setRoundingMode(FE_UPWARD);
    lowBdCoveringUp += other.lowBdCoveringUp;
    w1SumUP += other.w1SumUP;
    w2SumUP += other.w2SumUP;
    givenUp += other.givenUp;
    totalDensityAccountedUP += other.totalDensityAccountedUP;
    
    // Accumulate counting statistics
    count += other.count;
    givenUpcount += other.givenUpcount;
}

// ============================================================================
// CoveringMathFramework Implementation
// ============================================================================

void CoveringMathFramework::initializeCovering(int zpow, int numPrimes, int prec) {
    // Initialize base framework
    initialize(zpow, numPrimes, prec, "Cov");
    
    // Initialize Bell numbers cache
    initializeBellCache();
}

void CoveringMathFramework::initializeBellCache() {
    bellCache_ = CachedBellNums;
    ProgressUtils::logProgress("Initialized Bell numbers cache with " + 
                              std::to_string(bellCache_.size()) + " entries");
}

long long CoveringMathFramework::getBellNumber(int j, int k) const {
    // Handle special cases for small parameters
    if (j == 1) return k;
    if (j == 2) return -k * (k - 1);
    
    // Check bounds and retrieve from cache
    if (k > 0 && k <= static_cast<int>(bellCache_.size()) && 
        j > 0 && j <= static_cast<int>(bellCache_[k-1].size())) {
        return bellCache_[k-1][j-1];
    }
    
    // Comprehensive error message for out-of-bounds access
    std::ostringstream error_msg;
    error_msg << "FATAL ERROR: Bell number B(" << j << ", " << k << ") is out of cached range.\n"
              << "Cache dimensions: " << bellCache_.size() << " × " 
              << (bellCache_.empty() ? 0 : bellCache_[0].size()) << "\n"
              << "This indicates the computation has exceeded the precomputed Bell number table.\n"
              << "The computation cannot continue reliably without this value.\n"
              << "Solutions:\n"
              << "1. Reduce the number of primes parameter (numPrimes)\n"
              << "2. Reduce the search depth parameter (zpow)\n"
              << "3. Extend the Bell numbers cache in BellNums.h";
    
    throw std::runtime_error(error_msg.str());
}

// ============================================================================
// TraverseCoveringNumbers Implementation
// ============================================================================

inline void TraverseCoveringNumbers::newCoveringNumber(long double invSmoothUp, 
                                                      long double invSmoothDown, 
                                                      int k, 
                                                      SmoothNumbersCovVars& vars) const {
    // Process newly discovered covering number with interval arithmetic
    MathUtils::setRoundingMode(FE_DOWNWARD);
    vars.lowBdCovering += framework_.getDensityBounds(k - 1).first * invSmoothDown / framework_.getPrime(k);
    
    MathUtils::setRoundingMode(FE_UPWARD);
    const long double newup = framework_.getDensityBounds(k - 1).second * invSmoothUp / framework_.getPrime(k);
    vars.totalDensityAccountedUP += newup;
    vars.lowBdCoveringUp += newup;
}

inline void TraverseCoveringNumbers::newCoveringInterval(long double invSmoothUp, 
                                                        long double invSmoothDown, 
                                                        int initk, int k, 
                                                        SmoothNumbersCovVars& vars) const {
    // Process interval of covering numbers with rigorous bounds
    MathUtils::setRoundingMode(FE_DOWNWARD);
    vars.lowBdCovering += (framework_.getDensityBounds(initk - 1).first - 
                          framework_.getDensityBounds(k - 1).second) * invSmoothDown;
    
    MathUtils::setRoundingMode(FE_UPWARD);
    const long double newup = (framework_.getDensityBounds(initk - 1).second - 
                              framework_.getDensityBounds(k - 1).first) * invSmoothUp;
    vars.totalDensityAccountedUP += newup;
    vars.lowBdCoveringUp += newup;
}

long double TraverseCoveringNumbers::computeCoveringIndexDirect(const mpz_class& smoothNum, 
                                                               int sunDivCount, 
                                                               const mpz_class& ell, 
                                                               const std::vector<Divisor>& divisors) const {
    /**
     * Compute covering index using Bell number formula:
     * Formula: 1 + (ℓ-1)/ℓ + (1/ℓ) * Σ_{d|b, d>1} B(ω(d), τ(ℓ))/d
     * 
     * Note: Uses long double arithmetic for computational efficiency.
     * This should be sufficiently accurate for values that fit in 63-bit range,
     * but could be enhanced with arbitrary precision if needed.
     */
    long double ell_val = ell.get_d();
    long double covering_index = (2.0L * ell_val - 1.0L);  // Will be divided by ell
    
    for (const auto& div : divisors) {
        if (div.d > 1) {
            try {
                long long bell_val = framework_.getBellNumber(div.omega, sunDivCount);
                long double contribution = static_cast<long double>(bell_val) / div.d.get_d();
                covering_index += contribution;
            } catch (const std::runtime_error& e) {
                // Re-throw with additional context
                throw std::runtime_error("Error computing covering index for smooth number " + 
                                        smoothNum.get_str() + ": " + e.what());
            }
        }
    }
    
    return covering_index / ell_val;
}

inline long double TraverseCoveringNumbers::computeShiftedCoveringIndexSum(int sunDivCount,  
                                                                          const std::vector<Divisor>& divisors) const {
    /**
     * Compute shifted covering index sum for prime power iterations.
     * This precomputes a portion of the covering index that can be reused
     * across multiple prime power iterations.
     */
    long double covering_index = 0.0L; 
    
    for (const auto& div : divisors) {
        try {
            long long bell_val = framework_.getBellNumber(div.omega + 1, sunDivCount);
            covering_index += static_cast<long double>(bell_val) / div.d.get_d();
        } catch (const std::runtime_error& e) {
            throw std::runtime_error("Error computing shifted covering index: " + string(e.what()));
        }
    }
    
    return covering_index;
}

int TraverseCoveringNumbers::handleAutomaticCovering(int k, int sunDivCount, 
                                                    long double invSmoothUp, 
                                                    long double invSmoothDown, 
                                                    SmoothNumbersCovVars& vars, 
                                                    bool isDistrib) const {
    // Covering-specific: Handle automatic covering via small prime condition
    while (k <= framework_.numPrimes && framework_.getPrime(k) <= sunDivCount) {
        if (!isDistrib || vars.taskIndex == 0) {
            newCoveringNumber(invSmoothUp, invSmoothDown, k, vars);
        }
        k++;
        if (k == framework_.numPrimes + 1) {
            if (!isDistrib || vars.taskIndex == 0) {
                reachPrimeBoundary(invSmoothUp, COVERING_THRESHOLD, vars);  // Use covering threshold
            }
            return k;
        }
    }
    return k;
}

// Main computation control methods
int TraverseCoveringNumbers::determineTaskCount(SmoothNumbersCovVars& final) {
    ProgressUtils::logProgress("Determining parallel task decomposition");
    
    SmoothNumbersCovVars vars0;
    vars0.taskIndex = 0;
    std::vector<Divisor> initialDivisors = {{mpz_class(1), 0}};
    mpz_class initialEll(1);
    smoothNumbersCovDistrib(mpz_class(1), 1.0L, 1.0L, 1, 1, initialEll, initialDivisors, vars0);
    
    final.mergeWith(vars0);
    
    const int numTasks = vars0.givenUpcount;
    ProgressUtils::logProgress("Task decomposition complete. Tasks: " + std::to_string(numTasks));
    return numTasks;
}

void TraverseCoveringNumbers::loadExistingResults(const string& infilename, const string& outfilename,
                                                  std::vector<bool>& processedTasks, SmoothNumbersCovVars& final,
                                                  std::ofstream& fout) {
    if (!FileUtils::fileExists(infilename)) {
        return;
    }
    
    ProgressUtils::logProgress("Loading existing results from " + infilename);
    std::ifstream fin(infilename, std::ios::binary);
    SmoothNumbersCovVars ivars;
    
    int loadedCount = 0;
    while (fin.read(reinterpret_cast<char*>(&ivars), sizeof(ivars))) {
        if (fin.eof()) break;
        
        if (ivars.taskIndex > 0 && ivars.taskIndex <= static_cast<int>(processedTasks.size()) && 
            !processedTasks[ivars.taskIndex - 1]) {
            fout.write(reinterpret_cast<const char*>(&ivars), sizeof(ivars));
            fout.flush();
            final.mergeWith(ivars);
            processedTasks[ivars.taskIndex - 1] = true;
            loadedCount++;
        }
    }
    
    fin.close();
    ProgressUtils::logProgress("Loaded " + std::to_string(loadedCount) + " existing results");
}

void TraverseCoveringNumbers::computeRemainingTasks(const int numTasks, const std::vector<bool>& processedTasks,
                                                    SmoothNumbersCovVars& final, std::ofstream& fout) {
    ProgressUtils::logProgress("Computing remaining tasks in parallel");

    int remainingTasks = 0;
    for (int i = 0; i < numTasks; i++) {
        if (!processedTasks[i]) remainingTasks++;
    }
    
    ProgressUtils::logProgress("Remaining tasks: " + std::to_string(remainingTasks));

    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 1; i <= numTasks; i++) {
        if (processedTasks[i - 1]) continue;
        
        SmoothNumbersCovVars vars;
        vars.taskIndex = i;
        
        try {
            std::vector<Divisor> initialDivisors = {{mpz_class(1), 0}};
            mpz_class initialEll(1);
            smoothNumbersCovDistrib(mpz_class(1), 1.0L, 1.0L, 1, 1, initialEll, initialDivisors, vars);
            
            #pragma omp critical
            {
                ProgressUtils::logProgress("Completed task " + std::to_string(i));
                fout.write(reinterpret_cast<const char*>(&vars), sizeof(vars));
                fout.flush();
                final.mergeWith(vars);
            }
        } catch (const std::exception& e) {
            #pragma omp critical
            {
                cerr << "Error in task " << i << ": " << e.what() << endl;
            }
        }
    }
}

void TraverseCoveringNumbers::outputFinalResults(const SmoothNumbersCovVars& final) const {
    ProgressUtils::logProgress("Computation complete! Analyzing results...");
    
    // Validate total density accounting  
    if (!ValidationUtils::validateTotalDensity(final.totalDensityAccountedUP)) {
        cerr << "Warning: Total density not fully accounted for: " 
             << final.totalDensityAccountedUP << endl;
    }
    
    // Compute final bounds
    const long double upBdCovering = final.w1SumUP + final.w2SumUP + final.givenUp + final.lowBdCoveringUp;
    
    ProgressUtils::printReadableResults(
        "Covering",
        final.count,
        static_cast<double>(final.lowBdCovering),
        static_cast<double>(upBdCovering),
        static_cast<double>(final.w1SumUP),
        static_cast<double>(final.w2SumUP),
        static_cast<double>(final.totalDensityAccountedUP),
        static_cast<double>(final.givenUp),
        static_cast<double>(final.lowBdCoveringUp)
    );
}

void TraverseCoveringNumbers::calculate(const int zpow, const int numPrimes, const int prec) {
    try {
        // Initialize shared mathematical framework with covering-specific components
        framework_.initializeCovering(zpow, numPrimes, prec);
        
        // Initialize final results container
        SmoothNumbersCovVars final;
        
        // Determine parallel task decomposition
        const int numTasks = determineTaskCount(final);
        
        // Set up file I/O for persistence and recovery
        const string outfilename = FileUtils::generateDataFilename("SavedCovData", numPrimes, RLEN, zpow, "out.dat");
        const string infilename = FileUtils::generateDataFilename("SavedCovData", numPrimes, RLEN, zpow, "in.dat");
        
        std::ofstream fout(outfilename, std::ios::binary);
        if (!fout) {
            throw std::runtime_error("Cannot open output file: " + outfilename);
        }
        
        std::vector<bool> processedTasks(numTasks, false);
        
        // Load any existing results from previous runs
        loadExistingResults(infilename, outfilename, processedTasks, final, fout);
        
        // Compute remaining tasks in parallel
        computeRemainingTasks(numTasks, processedTasks, final, fout);
        
        fout.close();
        
        // Output final results
        outputFinalResults(final);
        
    } catch (const std::exception& e) {
        cerr << "\nFatal error in covering number calculation: " << e.what() << endl;
        throw;
    }
}

// ============================================================================
// Core Traversal Implementation
// ============================================================================

/**
 * Core covering number traversal algorithm (non-distributed version)
 * 
 * This function implements the main algorithm for computing covering number density bounds.
 * It recursively explores smooth numbers and computes their covering indices using Bell numbers
 * to determine if they contribute to the density calculation. In order to be a covering number, 
 * n must satisfy c'(n) >=2 which is what is tested for here.
 * 
 * @param smoothNum Current smooth number being processed
 * @param invSmoothUp Upper bound of 1/smoothNum for interval arithmetic
 * @param invSmoothDown Lower bound of 1/smoothNum for interval arithmetic
 * @param k Current prime index in the traversal
 * @param sunDivCount Current value of τ(ℓ) (divisor count parameter)
 * @param coveringIndex Current covering index (σ(n)/n value)
 * @param ell Current ℓ parameter for Bell number computation
 * @param divisors Current divisor structure for covering index calculation
 * @param vars Accumulation variables for results
 */
void TraverseCoveringNumbers::smoothNumbersCov(const mpz_class& smoothNum, 
                                              const long double invSmoothUp, 
                                              const long double invSmoothDown, 
                                              int k, int sunDivCount, 
                                              const long double coveringIndex,
                                              const mpz_class& ell,
                                              const std::vector<Divisor>& divisors, 
                                              SmoothNumbersCovVars& vars) {
    
    // Precompute shifted covering for prime power calculations
    long double shiftedCovering = computeShiftedCoveringIndexSum(sunDivCount, divisors);
    
    // Declare akReturn for reuse throughout the function
    long double akReturn;
    
    // Base case: exhausted prime set
    if (k == framework_.numPrimes + 1) {
        reachPrimeBoundary(invSmoothUp, coveringIndex, vars);
        return;
    }

    // Handle automatic covering detection
    k = handleAutomaticCovering(k, sunDivCount, invSmoothUp, invSmoothDown, vars, false);

    // Main loop over prime indices
    while (k <= framework_.numPrimes) {
        if (k == framework_.numPrimes + 1) {
            reachPrimeBoundary(invSmoothUp, coveringIndex, vars);
            return;
        }

        // Check for termination condition
        if (shouldTerminateAtPrimeLevel(coveringIndex, k, invSmoothDown)) {
            MathUtils::setRoundingMode(FE_DOWNWARD);
            auto [lambdaValues, lambdaRatios] = framework_.getLambdaBounds(k - 1);
            long double akReturn = framework_.akBound(AK_BOUND_TARGET, coveringIndex, lambdaValues, lambdaRatios);
            // Note: akBound sets rounding mode to FE_UPWARD before returning
            addToW1Sum(smoothNum, invSmoothUp, akReturn, k - 1, vars);
            return;
        }

        // Process current prime and its powers
        const int q = framework_.getPrime(k);
        int power = 0;
        mpz_class primePower(1);  

        mpz_class currentNum = smoothNum;
        long double currentInvUp = invSmoothUp;
        long double currentInvDown = invSmoothDown;
        long double newCovering = coveringIndex;
        int currentSunDivCount = sunDivCount;
        mpz_class currentEll = ell;
        std::vector<Divisor> currentDivisors = divisors;
        
        // Inner loop over prime powers
        while (true) {
            // Multiply by next prime power
            currentNum *= q;
            primePower *= q;
            power++;
            
            // Update covering index efficiently
            newCovering += shiftedCovering / (ell.get_d() * primePower.get_d());
            
            // Check if covering number found
            if (newCovering >= COVERING_THRESHOLD) {
                newCoveringNumber(currentInvUp, currentInvDown, k, vars);
                break;
            }
            
            // Update inverse values with proper rounding
            currentInvUp /= q;
            MathUtils::setRoundingMode(FE_DOWNWARD);
            currentInvDown /= q;
            MathUtils::setRoundingMode(FE_UPWARD);
            
            // Check bound for giving up on this prime
            if (shouldGiveUpOnPrimePower(newCovering, k, currentInvDown, false)) {
                // Add to W₂ sum
                long double finalCovering = coveringIndex + shiftedCovering / (ell.get_d() * (q - 1));
                MathUtils::setRoundingMode(FE_DOWNWARD);
                auto [lambdaValues, lambdaRatios] = framework_.getLambdaBounds(k);
                akReturn = framework_.akBound(AK_BOUND_TARGET, finalCovering, 
                                                           lambdaValues, lambdaRatios);
                // Note: akBound sets rounding mode to FE_UPWARD before returning
                
                incrementCount(1, currentNum, k, vars);
                auto densityBounds = framework_.getDensityBounds(k - 1);
                vars.totalDensityAccountedUP += densityBounds.second * currentInvUp;
                vars.w2SumUP += akReturn * densityBounds.second * currentInvUp;
                break;
            }

            // Covering-specific: Update parameters for special cases
            if (q == sunDivCount + 1) {
                currentSunDivCount = sunDivCount * (power + 1);
                currentEll = currentEll * q;
            } 
            
            if (q > sunDivCount + 1) {
                // Update divisors for covering computation
                for (const auto& div : divisors) {
                    currentDivisors.emplace_back(div.d * primePower, div.omega + 1);
                }
            }

            // Recursive call to continue traversal
            smoothNumbersCov(currentNum, currentInvUp, currentInvDown, k + 1, 
                            currentSunDivCount, newCovering, currentEll, currentDivisors, vars);
            
            vars.count++;
        }
        
        k++;
    }
}

/**
 * Distributed version of covering number traversal for parallel computation
 * 
 * This function implements task decomposition for parallel processing of covering number
 * density calculations. It explores the same algorithm as smoothNumbersCov but distributes
 * work across multiple tasks.
 * 
 * @param smoothNum Current smooth number being processed
 * @param invSmoothUp Upper bound of 1/smoothNum for interval arithmetic
 * @param invSmoothDown Lower bound of 1/smoothNum for interval arithmetic
 * @param k Current prime index in the traversal
 * @param sunDivCount Current value of τ(ℓ) (divisor count parameter)
 * @param ell Current ℓ parameter for Bell number computation
 * @param divisors Current divisor structure for covering index calculation
 * @param vars Accumulation variables for results (includes taskIndex for distribution)
 */
void TraverseCoveringNumbers::smoothNumbersCovDistrib(const mpz_class& smoothNum, 
                                                     const long double invSmoothUp, 
                                                     const long double invSmoothDown, 
                                                     int k, int sunDivCount,
                                                     const mpz_class& ell,
                                                     const std::vector<Divisor>& divisors, 
                                                     SmoothNumbersCovVars& vars) {
    // Compute covering index and precompute shifted covering
    long double coveringIndex = computeCoveringIndexDirect(smoothNum, sunDivCount, ell, divisors);
    long double shiftedCovering = computeShiftedCoveringIndexSum(sunDivCount, divisors);

    if (k == framework_.numPrimes + 1) {
        if (vars.taskIndex == 0) {
            reachPrimeBoundary(invSmoothUp, coveringIndex, vars);
        }
        return;
    }

    // Handle automatic covering detection
    k = handleAutomaticCovering(k, sunDivCount, invSmoothUp, invSmoothDown, vars, true);

    // Main distribution loop
    while (k <= framework_.numPrimes) {
        if (k == framework_.numPrimes + 1) {
            if (vars.taskIndex == 0) {
                reachPrimeBoundary(invSmoothUp, coveringIndex, vars);
            }
            return;
        }
        
        // Compute the AK bound for termination decision
        // This determines if we should continue processing or terminate
        long double akReturn;
        if (k < 1) {
            // For initial cases, use default bound of 1.0
            akReturn = 1.0L;
        } else {
            // Calculate the AK bound using lambda arrays and current covering index
            MathUtils::setRoundingMode(FE_DOWNWARD);
            auto [lambdaValues, lambdaRatios] = framework_.getLambdaBounds(k - 1);
            akReturn = framework_.akBound(AK_BOUND_TARGET, coveringIndex, lambdaValues, lambdaRatios);
            // Note: akBound sets rounding mode to FE_UPWARD before returning
        }

        if (k > 0 && shouldTerminateAtPrimeLevel(coveringIndex, k, invSmoothDown)) {
            if (vars.taskIndex == 0) {
                addToW1Sum(smoothNum, invSmoothUp, akReturn, k - 1, vars);
            }
            return;
        }
        
        if (k > 0) {
            // Process prime powers with task distribution
            const int q = framework_.getPrime(k);
            mpz_class currentNum = smoothNum;
            long double currentInvUp = invSmoothUp;
            long double currentInvDown = invSmoothDown;
            long double newCovering = coveringIndex;
            int currentSunDivCount = sunDivCount;
            mpz_class currentEll = ell;
            std::vector<Divisor> currentDivisors = divisors;
            int power = 0;
            mpz_class primePower(1);
            
            // Loop over powers of this prime
            while (true) {
                currentNum *= q;
                power++;
                primePower *= q;
                
                newCovering += shiftedCovering / (ell.get_d() * primePower.get_d());

                if (newCovering >= COVERING_THRESHOLD) {
                    if (vars.taskIndex == 0) newCoveringNumber(currentInvUp, currentInvDown, k, vars);
                    break;
                }

                currentInvUp /= q;
                MathUtils::setRoundingMode(FE_DOWNWARD);
                currentInvDown /= q;

                // Task distribution logic
                if (shouldGiveUpOnPrimePower(newCovering, k, currentInvDown, true)) {
                    vars.givenUpcount++;
                    if (vars.taskIndex == vars.givenUpcount) {
                        ProgressUtils::logProgress("Executing task " + std::to_string(vars.taskIndex) + 
                                                 " for number " + currentNum.get_str());
                        
                        // Execute full computation for this task
                        while (true) {
                            MathUtils::setRoundingMode(FE_DOWNWARD);
                            auto [lambdaValues, lambdaRatios] = framework_.getLambdaBounds(k);
                            akReturn = framework_.akBound(AK_BOUND_TARGET, newCovering, 
                                                                        lambdaValues, lambdaRatios);
                            // Note: akBound sets rounding mode to FE_UPWARD before returning
                            
                            if (akReturn * framework_.getDensityBounds(k).first * framework_.z * currentInvDown < k + 1) {
                                // Add to W₂ sum
                                long double finalCovering = coveringIndex + shiftedCovering / (ell.get_d() * (q - 1));
                                MathUtils::setRoundingMode(FE_DOWNWARD);
                                auto [finalLambdaValues, finalLambdaRatios] = framework_.getLambdaBounds(k);
                                akReturn = framework_.akBound(AK_BOUND_TARGET, finalCovering, 
                                                                           finalLambdaValues, finalLambdaRatios);
                                // Note: akBound sets rounding mode to FE_UPWARD before returning
                                
                                incrementCount(1, currentNum, k, vars);
                                auto densityBounds = framework_.getDensityBounds(k - 1);
                                vars.totalDensityAccountedUP += densityBounds.second * currentInvUp;
                                vars.w2SumUP += akReturn * densityBounds.second * currentInvUp;
                                break;
                            }
                            
                            // Update covering-specific parameters
                            if (q == sunDivCount + 1) {
                                currentSunDivCount = sunDivCount * (power + 1);
                                currentEll = ell * primePower;
                            }

                            if (q > sunDivCount + 1) {
                                for (const auto& div : divisors) {
                                    currentDivisors.emplace_back(div.d * primePower, div.omega + 1);
                                }
                            }

                            smoothNumbersCov(currentNum, currentInvUp, currentInvDown, k + 1, 
                                           currentSunDivCount, newCovering, currentEll, currentDivisors, vars);

                            // Move to next power
                            currentNum *= q;
                            primePower *= q;
                            newCovering += shiftedCovering / (ell.get_d() * primePower.get_d());
                            
                            if (newCovering >= COVERING_THRESHOLD) {
                                newCoveringNumber(currentInvUp, currentInvDown, k, vars);
                                break;
                            }

                            currentInvUp /= q;
                            MathUtils::setRoundingMode(FE_DOWNWARD);
                            currentInvDown /= q;
                            MathUtils::setRoundingMode(FE_UPWARD);
                            power++;

                            vars.count++;
                        }
                        return;
                    }
                    break;
                }
                
                // Update covering-specific parameters for distribution
                if (q == sunDivCount + 1) {
                    currentSunDivCount = sunDivCount * (power + 1);
                    currentEll = currentEll * q;
                }
                
                if (q > sunDivCount + 1) {
                    for (const auto& div : divisors) {
                        currentDivisors.emplace_back(div.d * primePower, div.omega + 1);
                    }
                }
                
                smoothNumbersCovDistrib(currentNum, currentInvUp, currentInvDown, k + 1, 
                                      currentSunDivCount, currentEll, currentDivisors, vars);
            }
        }
        
        k++;
    }
}

// ============================================================================
// Main Function
// ============================================================================

int main(int argc, char* argv[]) {
    cout << "Covering Number Density Calculator v19" << endl;
    cout << "====================================================" << endl;
    
    if (argc != 3 && argc != 4) {
        cout << "\nThis program computes rigorous upper and lower bounds for the density of covering numbers." << endl;
        cout << "\nUsage:" << endl;
        cout << "  " << argv[0] << " <zpow> <numPrimes> [precision]" << endl;
        cout << "\nParameters:" << endl;
        cout << "  zpow      - Search depth parameter (recommended: 10-20)" << endl;
        cout << "  numPrimes - Number of consecutive primes to consider (recommended: 50000-200000)" << endl;
        cout << "  precision - MPFR precision in bits (optional, default: 80, recommended: 80-120)" << endl;
        cout << "\nExamples:" << endl;
        cout << "  " << argv[0] << " 15 100000     # Uses default precision of 80 bits" << endl;
        cout << "  " << argv[0] << " 15 100000 120 # Uses custom precision of 120 bits" << endl;
        cout << "\nNote: This computation requires precomputed Bell numbers. The program will" << endl;
        cout << "      terminate with an error if the computation exceeds the cached range." << endl;
        return 0;
    }
    
    try {
        const int zpow = std::atoi(argv[1]);
        const int numPrimes = std::atoi(argv[2]);
        const int prec = (argc == 4) ? std::atoi(argv[3]) : 80;

        // Set initial rounding mode for interval arithmetic
        MathUtils::setRoundingMode(FE_UPWARD);
        
        TraverseCoveringNumbers calculator;
        calculator.calculate(zpow, numPrimes, prec);
        
        return 0;
        
    } catch (const std::runtime_error& e) {
        cout << "\n" << std::string(60, '=') << endl;
        cout << "COMPUTATION TERMINATED DUE TO ERROR:" << endl;
        cout << std::string(60, '=') << endl;
        cout << e.what() << endl;
        cout << std::string(60, '=') << endl;
        cout << "\nPossible solutions:" << endl;
        cout << "1. Reduce the number of primes parameter (numPrimes)" << endl;
        cout << "2. Reduce the search depth parameter (zpow)" << endl;
        cout << "3. Extend the Bell numbers cache in BellNums.h" << endl;
        cout << std::string(60, '=') << endl;
        return 1;
    } catch (const std::exception& e) {
        cerr << "\nError: " << e.what() << endl;
        return 1;
    } catch (...) {
        cerr << "\nUnknown error occurred" << endl;
        return 1;
    }
}
