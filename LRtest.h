#ifndef LR_UNENC
#define LR_UNENC

#include <NTL/ZZ.h>
#include "Database.h"
#include "../src/CZZ.h"   // -> for comparison

using namespace std;
using namespace NTL;

ZZ scaleup(double value, const long& l);
double scaledown(double value, const long& l);

// structure for LRparam
typedef struct LRpar{
    
    int max_iter;
    
    int dim1;
    int n_training;
    int logn;     // log2(n_training)
    
    double coeff[10];  // double coefficients
    ZZ evalcoeff[10];  // transformed coefficients for HE
    
    int polyscale;
    int log2polyscale;
    
    int polydeg;
    
    long logp;
    long nslots;
    
}LRpar;


void ReadLRparams(LRpar& LRparams, int max_iter, dMat zTrain,  int polydeg, long logp);


//---------------------------------------------------------------------
//! @ Plaintext logistic regression
void LR_sigmoid(dVec& theta, dMat zTrain, LRpar& LRparams);

void LR_poly(dVec& theta, dMat zTrain, LRpar& LRparams);


//---------------------------------------------------------------------
double getAUC(dVec theta, dMat zTest);

double getMSE(dVec theta1, dVec theta2);

double getNMSE(dVec theta1, dVec theta2);

#endif
