#ifndef SCHEME_PARAMS_H_
#define SCHEME_PARAMS_H_

#include <NTL/ZZ.h>

#include "KsiPows.h"


using namespace std;
using namespace NTL;

class Params {
public:
    long M;
    long N;
    long Nh;
    long logN;
    long logNh;
 
    long logp;    // base modulus 
    long logq;    // the largest modulus in the scheme
    long logP;    // modulus for key-switching
    long logPq;
    
    double sigma;
    double rho;
    long h;
    
    ZZ p;
    ZZ q;
    ZZ Pq;
    
    //ZZ* qi;
    //ZZ* Pqi;
    
    /**
     * rotation group for rotating messages withing slots
     */
    long** rotGroup;
    long** rotGroupInv;
    
    //-----------------------------------------
    
    Params(long logN,  long logp, long logq, double sigma = 3.2, double rho = 0.5, long h = 64);
};

#endif /* SCHEME_PARAMS_H_ */
