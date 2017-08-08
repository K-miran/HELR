#ifndef DATABASE_H
#define DATABASE_H

#include <vector>
#include <NTL/ZZ.h>

using namespace NTL;
using namespace std;

typedef vector<double> dVec;
typedef vector< vector<double> > dMat;

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
    
    long Slots;
    
}LRpar;



//!@ return: number of samples
int readData(dMat& Z, char* filename);

//---------------------------------------------------------------------
//!@ Choose the learning data and test data

void RandomSamplingData(dMat& zTrain, dMat& Ztest, dMat& Z);

void SamplingData(dMat& zTrain, dMat& Ztest, dMat& Z);

void ReadLRparams(LRpar& LRparams, int max_iter, dMat zTrain, double* coeff, ZZ* evalcoeff,  int polyscale, int polydeg, long Slots);

#endif
