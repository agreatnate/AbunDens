/**
 * DensityMathFramework.cpp - Implementation of shared mathematical framework
 */

#include "DensityMathFramework.h"

DensityMathFramework::~DensityMathFramework() {
    if (lambdas_) {
        MemoryUtils::deallocate2DArray(lambdas_, numPrimes + 1);
        MemoryUtils::deallocate2DArray(lambdasRatio_, numPrimes + 1);
    }
}

void DensityMathFramework::initialize(int zpow, int numPrimes, int prec, const string& problemType) {
    if (!ValidationUtils::validateParameters(zpow, numPrimes, prec)) {
        throw std::invalid_argument("Invalid computation parameters");
    }
    
    ProgressUtils::logProgress("Initializing mathematical framework for " + problemType + " numbers");
    
    initializePrimes(numPrimes, prec);
    rbd = initializeLambdas(numPrimes, prec);
    z = std::pow(2.0L, zpow);
    zz = zpow * zpow * zpow;  // For distributed computation
    this->prec = prec;
    
    ProgressUtils::logProgress("Framework initialization complete. Search depth: " + 
                              std::to_string(z));
}

void DensityMathFramework::initializePrimes(const int numPrimes, const int prec) {
    ProgressUtils::logProgress("Initializing prime sequence and density bounds");
    this->numPrimes = numPrimes;
    
    // Allocate arrays for prime data
    ourPrimes_ = std::make_unique<int[]>(numPrimes + 1);
    ourDensitiesDup_ = std::make_unique<long double[]>(numPrimes + 1);
    ourDensitiesDdn_ = std::make_unique<long double[]>(numPrimes + 1);
    
    // Initialize with first "prime" as 1 for algorithmic convenience
    long long p = 1;
    ourPrimes_[0] = 1;
    mpz_class curP(1);
    mpq_class curDensity(1);  // Running product of (p-1)/p
    ourDensitiesDup_[0] = 1.0L;
    ourDensitiesDdn_[0] = 1.0L;
    
    // High-precision temporary for density computation
    mpfr_t tempDensity;
    mpfr_init2(tempDensity, prec);
    
    // Generate consecutive primes and compute density bounds
    for (int count = 1; count <= numPrimes; count++) {
        if (count % 10000 == 0) {
            ProgressUtils::logProgress("Computing prime " + std::to_string(count));
        }
        
        p = PrimeUtils::nextPrime(p);
        if (!ValidationUtils::checkOverflow(p, static_cast<long long>(INT_MAX), 
                                           "prime generation")) {
            throw std::runtime_error("Prime sequence exceeds maximum size");
        }
        
        mpz_set_si(curP.get_mpz_t(), p);
        ourPrimes_[count] = static_cast<int>(p);
        
        // Update running density product: ∏(p-1)/p
        curDensity *= mpq_class(curP - 1, curP);
        
        // Compute rigorous interval bounds for density
        mpfr_set_q(tempDensity, curDensity.get_mpq_t(), MPFR_RNDD);
        ourDensitiesDdn_[count] = mpfr_get_ld(tempDensity, MPFR_RNDD);
        mpfr_set_q(tempDensity, curDensity.get_mpq_t(), MPFR_RNDU);
        ourDensitiesDup_[count] = mpfr_get_ld(tempDensity, MPFR_RNDU);
    }
    
    maxPrime = static_cast<int>(p);
    
    // Build inverse prime counting function π⁻¹(x)
    piX_ = std::make_unique<int[]>(maxPrime + 1);
    int i = maxPrime;
    for (int count = numPrimes; count >= 0; count--) {
        while (i >= ourPrimes_[count]) {
            piX_[i] = count;
            i--;
        }
    }
    
    if (i >= 0) piX_[0] = 0;
    
    mpfr_clear(tempDensity);
    ProgressUtils::logProgress("Prime initialization complete. Largest prime: " + 
                              std::to_string(maxPrime));
}

void DensityMathFramework::eulerMultiplicandBound(const mpq_class& invP, 
                                                 const mpz_class& denom, 
                                                 mpfr_t* Lambda, 
                                                 const mpz_class& curP, 
                                                 const unsigned long long r, 
                                                 const int rpow, 
                                                 const int prec) const {
    // Thread-safe static variables for temporary computation
    static thread_local bool initialized = false;
    static thread_local mpfr_t tmp, tmp2, tmp3;
    static thread_local int currentPrec = 0;
    
    if (!initialized || currentPrec != prec) {
        if (initialized) {
            mpfr_clear(tmp);
            mpfr_clear(tmp2); 
            mpfr_clear(tmp3);
        }
        mpfr_init2(tmp, prec);
        mpfr_init2(tmp2, prec);
        mpfr_init2(tmp3, prec);
        initialized = true;
        currentPrec = prec;
    }
    
    // Compute Euler multiplicand bound using high-precision arithmetic
    // Formula involves geometric series and power computations
    mpfr_set_ui(tmp, 1, MPFR_RNDU);
    mpfr_add_q(tmp, tmp, invP.get_mpq_t(), MPFR_RNDU);
    mpfr_add_q(tmp3, tmp, mpq_class(invP * invP).get_mpq_t(), MPFR_RNDU);
    mpfr_pow_si(tmp, tmp, static_cast<long>(r), MPFR_RNDU);
    mpfr_pow_si(tmp3, tmp3, static_cast<long>(r), MPFR_RNDU);

    mpfr_mul_z(tmp, tmp, mpz_class(curP - 1).get_mpz_t(), MPFR_RNDU);
    mpfr_add(tmp, tmp, tmp3, MPFR_RNDU);
    mpfr_sub_z(tmp, tmp, curP.get_mpz_t(), MPFR_RNDU);
    mpfr_mul_q(tmp, tmp, invP.get_mpq_t(), MPFR_RNDU);
    mpfr_mul_q(tmp, tmp, invP.get_mpq_t(), MPFR_RNDU);
    mpfr_add_si(tmp, tmp, 1, MPFR_RNDU);

    mpfr_set_q(tmp2, mpq_class(curP, curP - 1).get_mpq_t(), MPFR_RNDU);
    mpfr_pow_ui(tmp2, tmp2, r - 1, MPFR_RNDU);
    mpfr_mul_ui(tmp2, tmp2, r, MPFR_RNDU);
    mpfr_div_z(tmp2, tmp2, denom.get_mpz_t(), MPFR_RNDU);
    mpfr_add(tmp2, tmp, tmp2, MPFR_RNDU);
    mpfr_mul(Lambda[rpow], Lambda[rpow], tmp2, MPFR_RNDU);
}

void DensityMathFramework::delegliseBound(const mpz_class& curP, mpfr_t* Lambda, const int prec) const {
    // Deleglise bound computation for prime curP
    const mpz_class denom = (curP * curP - 1) * (curP * curP * curP * curP);
    const mpq_class invP(1, curP);

    unsigned long long r = 1;
    for (int rpow = 0; rpow < RLEN; rpow++) {
        eulerMultiplicandBound(invP, denom, Lambda, curP, r, rpow, prec);
        r *= 2;  // Powers of 2 for r parameter
    }
}


int DensityMathFramework::initializeLambdas(const int numPrimes, const int prec) {
    ProgressUtils::logProgress("Initializing lambda bound arrays");
    
    // Allocate 2D arrays for lambda values and ratios
    lambdas_ = MemoryUtils::allocate2DArray(numPrimes + 1, RLEN);
    lambdasRatio_ = MemoryUtils::allocate2DArray(numPrimes + 1, RLEN);

    const string filename = "LamdasStore-" + std::to_string(numPrimes) + '-' + std::to_string(prec) + ".dat";
    
    // Try loading from existing file first
    if (FileUtils::fileExists(filename)) {
        ProgressUtils::logProgress("Loading lambda values from " + filename);
        std::ifstream fin(filename, std::ios::binary);
        
        for (int i = numPrimes; i >= 0; i--) {
            fin.read(reinterpret_cast<char*>(lambdas_[i]), RLEN * sizeof(long double));
            fin.read(reinterpret_cast<char*>(lambdasRatio_[i]), RLEN * sizeof(long double));
            
            if (fin.fail()) {
                throw std::runtime_error("Error reading lambda data from file: " + filename);
            }
        }
        
        fin.close();
        ProgressUtils::logProgress("Lambda values loaded successfully");
        return RLEN - 1;
    }

    // Compute lambda values from scratch
    ProgressUtils::logProgress("Computing lambda values (this may take a while)");
    std::ofstream fout(filename, std::ios::binary);
    if (!fout) {
        throw std::runtime_error("Cannot create lambda file: " + filename);
    }

    // Initialize high-precision lambda array
    std::unique_ptr<mpfr_t[]> Lambda = std::make_unique<mpfr_t[]>(RLEN);
    
    for (int i = 0; i < RLEN; i++) {
        mpfr_init2(Lambda[i], prec);
        // Initialize with exponential growth pattern
        mpfr_set_q(Lambda[i], mpq_class(7 * (1 << i), 10000000000).get_mpq_t(), MPFR_RNDU);
        mpfr_exp(Lambda[i], Lambda[i], MPFR_RNDU);
    }

    // Compute bounds for large primes beyond our working set
    mpz_class curP = ourPrimes_[numPrimes];
    while (true) {
        mpz_nextprime(curP.get_mpz_t(), curP.get_mpz_t());
        if (curP > 100000000) break;  // Practical upper limit
        delegliseBound(curP, Lambda.get(), prec);
    }
    
    // Compute and store lambda values for each prime in our working set
    for (int i = numPrimes; i >= 0; i--) {
        if (i % 1000 == 0) {
            ProgressUtils::logProgress("Computing lambda bounds for prime index " + std::to_string(i));
        }
        
        curP = ourPrimes_[i];
        
        for (int rpow = 0; rpow < RLEN; rpow++) {
            lambdas_[i][rpow] = mpfr_get_ld(Lambda[rpow], MPFR_RNDU) - 1.0L;
            lambdasRatio_[i][rpow] = (rpow == 0) ? 
                lambdas_[i][rpow] - 1.0L : 
                lambdas_[i][rpow] / lambdas_[i][rpow - 1] - 1.0L;
        }
        
        // Write to file for future use
        fout.write(reinterpret_cast<const char*>(lambdas_[i]), RLEN * sizeof(long double));
        fout.write(reinterpret_cast<const char*>(lambdasRatio_[i]), RLEN * sizeof(long double));
        
        if (fout.fail()) {
            throw std::runtime_error("Error writing lambda data to file: " + filename);
        }
        
        // Update bounds for next iteration
        delegliseBound(curP, Lambda.get(), prec);
    }
    
    fout.close();
    
    // Clean up high-precision temporaries
    for (int i = 0; i < RLEN; i++) {
        mpfr_clear(Lambda[i]);
    }
    
    ProgressUtils::logProgress("Lambda computation complete and saved to " + filename);
    return RLEN - 1;
}
