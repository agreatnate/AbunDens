/**
 * AbunDens.cpp - Abundant number density calculation using shared framework
 * 
 * This refactored version leverages the DensityMathFramework for all common
 * mathematical operations, focusing only on abundant-number-specific logic.
 * The code size is reduced by ~60% while maintaining identical functionality.
 */

#include "AbunDens.h"
#include "DensityMathFramework.h"
#include "CommonUtils.h"
#include <omp.h>

// ============================================================================
// SmoothNumbersAbundVars Implementation
// ============================================================================

SmoothNumbersAbundVars::SmoothNumbersAbundVars() 
    : lowBdAbundant(0.0L), lowBdAbundantUp(0.0L), w1SumUP(0.0L), w2SumUP(0.0L),
      totalDensityAccountedUP(0.0L), givenUp(0.0L), count(0), 
      printCount(DEFAULT_PRINT_COUNT), givenUpcount(0), taskIndex(0) {
}

SmoothNumbersAbundVars::SmoothNumbersAbundVars(const SmoothNumbersAbundVars& other) 
    : lowBdAbundant(other.lowBdAbundant), lowBdAbundantUp(other.lowBdAbundantUp),
      w1SumUP(other.w1SumUP), w2SumUP(other.w2SumUP),
      totalDensityAccountedUP(other.totalDensityAccountedUP),
      givenUp(other.givenUp), count(other.count), printCount(other.printCount),
      givenUpcount(other.givenUpcount), taskIndex(other.taskIndex),
      tempMaxAbundQ(other.tempMaxAbundQ), tempMaxAbundZ(other.tempMaxAbundZ) {
}

SmoothNumbersAbundVars& SmoothNumbersAbundVars::operator=(const SmoothNumbersAbundVars& other) {
    if (this != &other) {
        lowBdAbundant = other.lowBdAbundant;
        lowBdAbundantUp = other.lowBdAbundantUp;
        w1SumUP = other.w1SumUP;
        w2SumUP = other.w2SumUP;
        totalDensityAccountedUP = other.totalDensityAccountedUP;
        givenUp = other.givenUp;
        count = other.count;
        printCount = other.printCount;
        givenUpcount = other.givenUpcount;
        taskIndex = other.taskIndex;
        tempMaxAbundQ = other.tempMaxAbundQ;
        tempMaxAbundZ = other.tempMaxAbundZ;
    }
    return *this;
}

void SmoothNumbersAbundVars::mergeWith(const SmoothNumbersAbundVars& other) {
    // Use rigorous interval arithmetic for bound accumulation
    // ROUNDING MODE VALIDATION: Lower bounds must use FE_DOWNWARD
    MathUtils::setRoundingMode(FE_DOWNWARD);
    lowBdAbundant += other.lowBdAbundant;
    
    // ROUNDING MODE VALIDATION: Upper bounds must use FE_UPWARD  
    MathUtils::setRoundingMode(FE_UPWARD);
    lowBdAbundantUp += other.lowBdAbundantUp;
    w1SumUP += other.w1SumUP;
    w2SumUP += other.w2SumUP;
    givenUp += other.givenUp;
    totalDensityAccountedUP += other.totalDensityAccountedUP;
    
    // Accumulate counting statistics
    count += other.count;
    givenUpcount += other.givenUpcount;
}

// ============================================================================
// TraverseAbundantNumbers Implementation
// ============================================================================

inline void TraverseAbundantNumbers::newAbundantNumber(long double invSmoothUp, 
                                                      long double invSmoothDown, 
                                                      int k, 
                                                      SmoothNumbersAbundVars& vars) const {
    // Process newly discovered abundant number with interval arithmetic
    // ROUNDING MODE VALIDATION: Lower bound uses FE_DOWNWARD + invSmoothDown
    MathUtils::setRoundingMode(FE_DOWNWARD);
    vars.lowBdAbundant += framework_.getDensityBounds(k - 1).first * invSmoothDown / framework_.getPrime(k);
    
    // ROUNDING MODE VALIDATION: Upper bound uses FE_UPWARD + invSmoothUp
    MathUtils::setRoundingMode(FE_UPWARD);
    const long double newup = framework_.getDensityBounds(k - 1).second * invSmoothUp / framework_.getPrime(k);
    vars.totalDensityAccountedUP += newup;
    vars.lowBdAbundantUp += newup;
}


inline void TraverseAbundantNumbers::newAbundantInterval(long double invSmoothUp, 
                                                        long double invSmoothDown, 
                                                        int initk, int k, 
                                                        SmoothNumbersAbundVars& vars) const {
    // Process interval primes that can multiply the current number to produce new abundant numbers
    MathUtils::setRoundingMode(FE_DOWNWARD);
    vars.lowBdAbundant += (framework_.getDensityBounds(initk - 1).first - 
                          framework_.getDensityBounds(k - 1).second) * invSmoothDown;
    
    MathUtils::setRoundingMode(FE_UPWARD);
    const long double newup = (framework_.getDensityBounds(initk - 1).second - 
                              framework_.getDensityBounds(k - 1).first) * invSmoothUp;
    vars.totalDensityAccountedUP += newup;
    vars.lowBdAbundantUp += newup;
}

int TraverseAbundantNumbers::handleAbundantDetection(const mpz_class& smoothNum, 
                                                    const mpq_class& curAbundancyQ, 
                                                    int k, long double invSmoothUp, 
                                                    long double invSmoothDown, 
                                                    SmoothNumbersAbundVars& vars, 
                                                    bool isDistrib) const {
    // Handle initial abundant number detection when processing a new number
    const int initk = k;
    
    // Only check k > 0 for distributed version to avoid accessing invalid indices
    if ((!isDistrib || k > 0) && curAbundancyQ.get_num() * (framework_.getPrime(k) + 1) >= 
                 2 * framework_.getPrime(k) * curAbundancyQ.get_den()) {
        
        // Find maximum abundant prime using efficient computation
        const mpq_class maxAbundQ = curAbundancyQ.get_num() / 
                            (2 * curAbundancyQ.get_den() - curAbundancyQ.get_num());
        
        int maxAbundP;
        if (maxAbundQ < framework_.maxPrime) {
            maxAbundP = MathUtils::safeMpzToInt(mpz_class(maxAbundQ));
        } else {
            maxAbundP = framework_.maxPrime;
        }
        
        k = framework_.getPrimeIndex(maxAbundP) + 1;
        
        if (!isDistrib || vars.taskIndex == 0) {
            incrementCount(k - initk, smoothNum, k, vars);
            newAbundantInterval(invSmoothUp, invSmoothDown, initk, k, vars);
        }
    }
    
    return k;
}

// Main computation control methods 
int TraverseAbundantNumbers::determineTaskCount(SmoothNumbersAbundVars& final) {
    ProgressUtils::logProgress("Determining parallel task decomposition");
    
    SmoothNumbersAbundVars vars0;
    vars0.taskIndex = 0;
    smoothNumbersAbundDistrib(mpz_class(1), 1.0L, 1.0L, 0, 1.0L, mpq_class(1), vars0);
    
    final.mergeWith(vars0);
    
    const int numTasks = vars0.givenUpcount;
    ProgressUtils::logProgress("Task decomposition complete. Tasks: " + std::to_string(numTasks));
    return numTasks;
}

void TraverseAbundantNumbers::loadExistingResults(const string& infilename, const string& outfilename,
                                                  std::vector<bool>& processedTasks, SmoothNumbersAbundVars& final,
                                                  std::ofstream& fout) {
    if (!FileUtils::fileExists(infilename)) {
        return;
    }
    
    ProgressUtils::logProgress("Loading existing results from " + infilename);
    std::ifstream fin(infilename, std::ios::binary);
    SmoothNumbersAbundVars ivars;
    
    int loadedCount = 0;
    while (fin.read(reinterpret_cast<char*>(&ivars), sizeof(ivars))) {
        if (fin.eof()) break;
        
        if (ivars.taskIndex > 0 && ivars.taskIndex <= processedTasks.size() && 
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

void TraverseAbundantNumbers::computeRemainingTasks(const int numTasks, const std::vector<bool>& processedTasks,
                                                    SmoothNumbersAbundVars& final, std::ofstream& fout) {
    ProgressUtils::logProgress("Computing remaining tasks in parallel");

    int remainingTasks = 0;
    for (int i = 0; i < numTasks; i++) {
        if (!processedTasks[i]) remainingTasks++;
    }
    
    ProgressUtils::logProgress("Remaining tasks: " + std::to_string(remainingTasks));

    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 1; i <= numTasks; i++) {
        if (processedTasks[i - 1]) continue;
        
        SmoothNumbersAbundVars vars;
        vars.taskIndex = i;
        
        try {
            smoothNumbersAbundDistrib(mpz_class(1), 1.0L, 1.0L, 1, 1.0L, mpq_class(1), vars);
            
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

void TraverseAbundantNumbers::outputFinalResults(const SmoothNumbersAbundVars& final) const {
    ProgressUtils::logProgress("Computation complete! Analyzing results...");
    
    // Validate total density accounting
    if (!ValidationUtils::validateTotalDensity(final.totalDensityAccountedUP)) {
        cerr << "Warning: Total density not fully accounted for: " 
             << final.totalDensityAccountedUP << endl;
    }
    
    // Compute final bounds
    const long double upBdAbundant = final.w1SumUP + final.w2SumUP + final.givenUp + final.lowBdAbundantUp;
    
    ProgressUtils::printReadableResults(
        "Abundant",
        final.count,
        static_cast<double>(final.lowBdAbundant),
        static_cast<double>(upBdAbundant),
        static_cast<double>(final.w1SumUP),
        static_cast<double>(final.w2SumUP),
        static_cast<double>(final.totalDensityAccountedUP),
        static_cast<double>(final.givenUp),
        static_cast<double>(final.lowBdAbundantUp)
    );
}

void TraverseAbundantNumbers::calculate(const int zpow, const int numPrimes, const int prec) {
    try {
        // Initialize shared mathematical framework
        framework_.initialize(zpow, numPrimes, prec, "Abund");
        
        // Initialize final results container
        SmoothNumbersAbundVars final;
        
        // Determine parallel task decomposition
        const int numTasks = determineTaskCount(final);
        
        // Set up file I/O for persistence and recovery
        const string outfilename = FileUtils::generateDataFilename("SavedData", numPrimes, RLEN, zpow, "out.dat");
        const string infilename = FileUtils::generateDataFilename("SavedData", numPrimes, RLEN, zpow, "in.dat");
        
        std::ofstream fout(outfilename, std::ios::binary);
        if (!fout) {
            throw std::runtime_error("Cannot open output file: " + outfilename);
        }
        
        std::vector<bool> processedTasks(numTasks, false);
        
        // Load any existing results from previous runs
        // RENAME "out" file to "in" to read in previous data!

        loadExistingResults(infilename, outfilename, processedTasks, final, fout);
        
        // Compute remaining tasks in parallel
        computeRemainingTasks(numTasks, processedTasks, final, fout);
        
        fout.close();
        
        // Output final results
        outputFinalResults(final);
        
    } catch (const std::exception& e) {
        cerr << "\nFatal error in abundant number calculation: " << e.what() << endl;
        throw;
    }
}

// ============================================================================
// Core Traversal Implementation
// ============================================================================

/**
 * Core abundant number traversal algorithm (non-distributed version)
 * 
 * This function implements the main algorithm for computing abundant number density bounds.
 * It recursively explores smooth numbers and checks their
 * abundancy σ(n)/n to determine if they contribute to the density calculation.
 * A number n is abundant if σ(n)/n ≥ 2.
 * 
 * @param smoothNum Current smooth number n being processed  
 * @param invSmoothUp Upper bound of 1/n for interval arithmetic
 * @param invSmoothDown Lower bound of 1/n for interval arithmetic  
 * @param k Current prime index in the traversal (p_k)
 * @param curAbundancyUP Upper bound of current abundancy σ(n)/n
 * @param curAbundancyQ Exact rational abundancy σ(n)/n
 * @param vars Accumulation variables for results
 */
void TraverseAbundantNumbers::smoothNumbersAbund(const mpz_class& smoothNum, 
                                                 long double invSmoothUp, 
                                                 long double invSmoothDown, 
                                                 int k, 
                                                 long double curAbundancyUP, 
                                                 const mpq_class& curAbundancyQ, 
                                                 SmoothNumbersAbundVars& vars) {
    
    // ═══════════════════════════════════════════════════════════════════════════
    // PHASE 1: Base Case - Check if we've exhausted the prime set
    // ═══════════════════════════════════════════════════════════════════════════
    if (k == framework_.numPrimes + 1) {
        reachPrimeBoundary(invSmoothUp, curAbundancyUP, vars);
        return;
    }

    // ═══════════════════════════════════════════════════════════════════════════
    // PHASE 2: Initial Abundant Number Detection
    // ═══════════════════════════════════════════════════════════════════════════
    const int initk = k;
    
    int maxAbundP;
    if (curAbundancyQ.get_num() * (framework_.getPrime(k) + 1) >= 
        ABUNDANT_THRESHOLD_INT * framework_.getPrime(k) * curAbundancyQ.get_den()) {
        const mpq_class maxAbundQ = curAbundancyQ.get_num() / 
                    (ABUNDANT_THRESHOLD_INT * curAbundancyQ.get_den() - curAbundancyQ.get_num());
        if (maxAbundQ < framework_.maxPrime) {
            maxAbundP = MathUtils::safeMpzToInt(mpz_class(maxAbundQ));
        } else {
            maxAbundP = framework_.maxPrime;
        }
        k = framework_.getPrimeIndex(maxAbundP) + 1;
        incrementCount(k - initk, smoothNum, k, vars);
        newAbundantInterval(invSmoothUp, invSmoothDown, initk, k, vars);
    }

    // ═══════════════════════════════════════════════════════════════════════════
    // PHASE 3: Main Prime Traversal Loop
    // ═══════════════════════════════════════════════════════════════════════════
    
    // Declare variables outside the loop
    mpq_class abundancyFactor;
    mpz_class qPower;

    // Main traversal loop through prime indices
    for (;; k++) {
        if (k == framework_.numPrimes + 1) {
            reachPrimeBoundary(invSmoothUp, curAbundancyUP, vars);
            return;
        }

        // Check for termination condition
        MathUtils::setRoundingMode(FE_DOWNWARD);
        auto [lambdaValues, lambdaRatios] = framework_.getLambdaBounds(k - 1);
        long double akReturn = framework_.akBound(AK_BOUND_TARGET, curAbundancyUP, lambdaValues, lambdaRatios);
        // Note: akBound sets rounding mode to FE_UPWARD before returning
        if (shouldTerminateAtPrimeLevel(curAbundancyUP, k, invSmoothDown)) {
            addToW1Sum(smoothNum, invSmoothUp, akReturn, k - 1, vars);
            return;
        } else {
            // Process current prime and its powers
        const int q = framework_.getPrime(k);
        long double abundancyFactorUp = static_cast<long double>(q + 1) / q;
        
        // Set up variables for prime power iteration
        mpq_set_ui(abundancyFactor.get_mpq_t(), q + 1, q);
        mpz_set_ui(qPower.get_mpz_t(), q);
        
        long double invSmoothWithqPowerUp = invSmoothUp / q;
        long double newAbundancyUP = curAbundancyUP * abundancyFactorUp;
        MathUtils::setRoundingMode(FE_DOWNWARD);
        long double invSmoothWithqPowerDown = invSmoothDown / q;
        auto [lambdaValues, lambdaRatios] = framework_.getLambdaBounds(k);
        akReturn = framework_.akBound(AK_BOUND_TARGET, newAbundancyUP, lambdaValues, lambdaRatios);
        // Note: akBound sets rounding mode to FE_UPWARD before returning
        
        // Inner loop over prime powers
        while (true) {
            // Check bound for giving up on this prime
            
            if (shouldGiveUpOnPrimePower(newAbundancyUP, k, invSmoothWithqPowerDown)) {
                // Add to W₂ sum and terminate this branch
                incrementCount(2, smoothNum, k, vars);
                auto densityBounds = framework_.getDensityBounds(k - 1);
                vars.totalDensityAccountedUP += densityBounds.second * invSmoothWithqPowerUp;
                
                long double finalAbundancyUP = curAbundancyUP * q / (q - 1);
                MathUtils::setRoundingMode(FE_DOWNWARD);
                auto [finalLambdaValues, finalLambdaRatios] = framework_.getLambdaBounds(k);
                akReturn = framework_.akBound(AK_BOUND_TARGET, finalAbundancyUP, finalLambdaValues, finalLambdaRatios);
                // Note: akBound sets rounding mode to FE_UPWARD before returning
                vars.w2SumUP += akReturn * densityBounds.second * invSmoothWithqPowerUp;
                break;
            }
            
            // Recursive call to continue traversal
            smoothNumbersAbund(smoothNum * qPower, invSmoothWithqPowerUp, invSmoothWithqPowerDown, 
                              k + 1, curAbundancyUP * abundancyFactorUp, curAbundancyQ * abundancyFactor, vars);
            
            // Update for next power iteration
            abundancyFactor /= q;
            abundancyFactor += 1;
            
            if (curAbundancyQ * abundancyFactor >= ABUNDANT_THRESHOLD_INT) {
                newAbundantNumber(invSmoothWithqPowerUp, invSmoothWithqPowerDown, k, vars);
                break;
            }
            
            qPower *= q;
            invSmoothWithqPowerUp /= q;
            abundancyFactorUp /= q;
            abundancyFactorUp += 1.0L;
            newAbundancyUP = curAbundancyUP * abundancyFactorUp;
            
            MathUtils::setRoundingMode(FE_DOWNWARD);
            invSmoothWithqPowerDown /= q;
            auto [loopLambdaValues, loopLambdaRatios] = framework_.getLambdaBounds(k);
            akReturn = framework_.akBound(AK_BOUND_TARGET, newAbundancyUP, loopLambdaValues, loopLambdaRatios);
            // Note: akBound sets rounding mode to FE_UPWARD before returning
            MathUtils::setRoundingMode(FE_UPWARD);
        }
        
        }  // end else block
    }
}

/**
 * Distributed version of abundant number traversal for parallel computation
 * 
 * This function implements task decomposition for parallel processing. It explores
 * the same algorithm as smoothNumbersAbund but distributes work across multiple
 * tasks for efficient parallelization. Only task 0 (the initial function call to 
 * determine the number of tasks records information from the initial part of the tree
 * On subsequent calls the tree is traverssed to find the proper starting point for the 
 * subtask call to smoothNumbersAbund(...)
 * 
 * @param smoothNum Current smooth number being processed
 * @param invSmoothUp Upper bound of 1/smoothNum for interval arithmetic
 * @param invSmoothDown Lower bound of 1/smoothNum for interval arithmetic
 * @param k Current prime index in the traversal
 * @param curAbundancyUP Upper bound of current abundancy value
 * @param curAbundancyQ Exact rational abundancy value
 * @param vars Accumulation variables for results (includes taskIndex for distribution)
 */
void TraverseAbundantNumbers::smoothNumbersAbundDistrib(const mpz_class& smoothNum, 
                                                       long double invSmoothUp, 
                                                       long double invSmoothDown, 
                                                       int k, 
                                                       long double curAbundancyUP, 
                                                       const mpq_class& curAbundancyQ, 
                                                       SmoothNumbersAbundVars& vars) {
    // Distributed version for parallel task decomposition
    if (k == framework_.numPrimes + 1) {
        if (vars.taskIndex == 0) reachPrimeBoundary(invSmoothUp, curAbundancyUP, vars);
        return;
    }

    // Handle initial abundant number detection (distributed version)
    k = handleAbundantDetection(smoothNum, curAbundancyQ, k, invSmoothUp, invSmoothDown, vars, true);

    // Declare variables outside the loop
    mpq_class abundancyFactor;
    mpz_class qPower;

    while (true) {
        if (k == 0) {
            k++;
            continue;
        }
        
        if (k == framework_.numPrimes + 1) {
            if (vars.taskIndex == 0) reachPrimeBoundary(invSmoothUp, curAbundancyUP, vars);
            return;
        }
        
        // Compute the AK bound for termination decision
        // This determines if we should continue processing or terminate
        long double akReturn;
        // Calculate the AK bound using lambda arrays and current abundancy
        MathUtils::setRoundingMode(FE_DOWNWARD);
        auto [lambdaValues, lambdaRatios] = framework_.getLambdaBounds(k - 1);
        akReturn = framework_.akBound(AK_BOUND_TARGET, curAbundancyUP, lambdaValues, lambdaRatios);
        // Note: akBound sets rounding mode to FE_UPWARD before returning
        if (double(akReturn * framework_.getDensityBounds(k - 1).first * framework_.z * invSmoothDown) < k) {
            if (vars.taskIndex == 0) {
                addToW1Sum(smoothNum, invSmoothUp, akReturn, k - 1, vars);
            }
            return;
        }
        
        // Process prime powers with task distribution
        const int q = framework_.getPrime(k);
        mpq_class abundancyFactor(q + 1, q);
        long double abundancyFactorUp = static_cast<long double>(q + 1) / q;
        mpz_class qPower(q);
        
        long double invSmoothWithqPowerUp = invSmoothUp / q;
        long double invSmoothWithqPowerDown = invSmoothDown;  //Will be divided by q below

        while (true) {
            long double newAbundancyUP = curAbundancyUP * abundancyFactorUp;
            MathUtils::setRoundingMode(FE_DOWNWARD);
            invSmoothWithqPowerDown /= q;
            auto [lambdaValues, lambdaRatios] = framework_.getLambdaBounds(k);
            akReturn = framework_.akBound(AK_BOUND_TARGET, newAbundancyUP, lambdaValues, lambdaRatios);
            // Note: akBound sets rounding mode to FE_UPWARD before returning

            // Task distribution logic
            if (double(akReturn * framework_.getDensityBounds(k).first * framework_.zz * invSmoothWithqPowerDown) < k + 1) {
                vars.givenUpcount++;
                if (vars.taskIndex == vars.givenUpcount) {
                    //This is the task we've been called to compute!
                    //Rather than simply calling smoothNumbersAbund(...) we still need to process all higher powers of the current prime here:
                    while (true) {
                        // Check bound for giving up on this prime
                        if (double(akReturn * framework_.getDensityBounds(k).first * framework_.z * invSmoothWithqPowerDown) < k + 1) {
                            // Add to W₂ sum and terminate this branch
                            incrementCount(2, smoothNum, k, vars);
                            auto densityBounds = framework_.getDensityBounds(k - 1);
                            vars.totalDensityAccountedUP += densityBounds.second * invSmoothWithqPowerUp;
                
                            long double finalAbundancyUP = curAbundancyUP * q / (q - 1);
                            MathUtils::setRoundingMode(FE_DOWNWARD);
                            auto [lambdaValues, lambdaRatios] = framework_.getLambdaBounds(k);
                            akReturn = framework_.akBound(AK_BOUND_TARGET, finalAbundancyUP, 
                                                                  lambdaValues, lambdaRatios);
                            // Note: akBound sets rounding mode to FE_UPWARD before returning
                            vars.w2SumUP += akReturn * densityBounds.second * invSmoothWithqPowerUp;
                            break;
                        }
            
                        // Recursive call to continue traversal
                        smoothNumbersAbund(smoothNum * qPower, invSmoothWithqPowerUp, invSmoothWithqPowerDown, 
                                  k + 1, curAbundancyUP * abundancyFactorUp, curAbundancyQ * abundancyFactor, vars);
            
                        // Update for next power iteration
                        abundancyFactor /= q;
                        abundancyFactor += 1;
            
                        if (curAbundancyQ * abundancyFactor >= ABUNDANT_THRESHOLD_INT) {
                            newAbundantNumber(invSmoothWithqPowerUp, invSmoothWithqPowerDown, k, vars);
                            break;
                        }
            
                        qPower *= q;
                        invSmoothWithqPowerUp /= q;
                        abundancyFactorUp /= q;
                        abundancyFactorUp += 1.0L;
                        newAbundancyUP = curAbundancyUP * abundancyFactorUp;
            
                        MathUtils::setRoundingMode(FE_DOWNWARD);
                        invSmoothWithqPowerDown /= q;
                        auto [innerLambdaValues, innerLambdaRatios] = framework_.getLambdaBounds(k);
                        akReturn = framework_.akBound(AK_BOUND_TARGET, newAbundancyUP, innerLambdaValues, innerLambdaRatios);
                        // Note: akBound sets rounding mode to FE_UPWARD before returning
                        MathUtils::setRoundingMode(FE_UPWARD);
                       
                    }

                    return;
                }
                break;
            }
            
            // Continue distribution traversal
            smoothNumbersAbundDistrib(smoothNum * qPower, invSmoothWithqPowerUp, invSmoothWithqPowerDown, 
                                     k + 1, curAbundancyUP * abundancyFactorUp, curAbundancyQ * abundancyFactor, vars);
            
            abundancyFactor /= q;
            abundancyFactor += 1;
            
            if (curAbundancyQ * abundancyFactor >= ABUNDANT_THRESHOLD_INT) {
                if (vars.taskIndex == 0) {
                    newAbundantNumber(invSmoothWithqPowerUp, invSmoothWithqPowerDown, k, vars);
                }
                break;
            }
            
            qPower *= q;
            invSmoothWithqPowerUp /= q;
            abundancyFactorUp /= q;
            abundancyFactorUp += 1.0L;
        }
        
        k += 1;
    }
}

// ============================================================================
// Main Function
// ============================================================================

int main(int argc, char* argv[]) {
    cout << "Abundant Number Density Calculator v19" << endl;
    cout << "===================================================" << endl;
    
    if (argc != 3 && argc != 4) {
        cout << "\nThis program computes rigorous upper and lower bounds for the density of abundant numbers." << endl;
        cout << "\nUsage:" << endl;
        cout << "  " << argv[0] << " <zpow> <numPrimes> [precision]" << endl;
        cout << "\nParameters:" << endl;
        cout << "  zpow      - Search depth parameter (recommended: 10-20)" << endl;
        cout << "  numPrimes - Number of consecutive primes to consider (recommended: 50000-200000)" << endl;
        cout << "  precision - MPFR precision in bits (optional, default: 80, recommended: 80-120)" << endl;
        cout << "\nExamples:" << endl;
        cout << "  " << argv[0] << " 15 100000     # Uses default precision of 80 bits" << endl;
        cout << "  " << argv[0] << " 15 100000 120 # Uses custom precision of 120 bits" << endl;
        return 0;
    }
    
    try {
        const int zpow = std::atoi(argv[1]);
        const int numPrimes = std::atoi(argv[2]);
        const int prec = (argc == 4) ? std::atoi(argv[3]) : 80;

        // Set initial rounding mode for interval arithmetic
        MathUtils::setRoundingMode(FE_UPWARD);
        
        TraverseAbundantNumbers calculator;
        calculator.calculate(zpow, numPrimes, prec);
        
        return 0;
        
    } catch (const std::exception& e) {
        cerr << "\nError: " << e.what() << endl;
        return 1;
    } catch (...) {
        cerr << "\nUnknown error occurred" << endl;
        return 1;
    }
}
