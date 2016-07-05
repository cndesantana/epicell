#ifndef _RNG_H_
#define _RNG_H_

#include<random>

typedef std::mt19937                     NENG;    // Mersenne Twister
typedef std::uniform_int_distribution<> iDIST;   // Uniform Integer Distribution
typedef std::uniform_real_distribution<> fDIST;   // Uniform Real Distribution
typedef std::mt19937                     ENG;    // Mersenne Twister
typedef std::normal_distribution<> nDIST;   // Normal Distribution

class RNG {
private:
    static ENG eng;
    iDIST idist;
    fDIST fdist;
public:
    static void seed(int s) { eng.seed(s); }
    RNG(int imin, int imax) : idist(imin, imax) {}
    RNG(float fmin, float fmax) : fdist(fmin, fmax) {}
    int igenerate() { return idist(eng); }
    float fgenerate() { return fdist(eng); }
};

class NRNG {//random number generator from a Normal Distribution
private:
    static NENG neng;
    nDIST ndist;
public:
    static void seed(int s) { neng.seed(s); }
    NRNG(double mean, double sd) : ndist(mean, sd) {}
    double ngenerate() { return ndist(neng); }
};

// the one generator, stored as RNG::eng
ENG RNG::eng;
// the one generator, stored as NRNG::neng
NENG NRNG::neng;

#endif
