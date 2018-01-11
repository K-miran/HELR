#ifndef DATABASE_H
#define DATABASE_H

#include <vector>
#include <NTL/ZZ.h>

using namespace NTL;
using namespace std;

typedef vector<double> dVec;
typedef vector< vector<double> > dMat;




//!@ return: number of samples
int readData(dMat& Z, char* filename);

//---------------------------------------------------------------------
//!@ Choose the learning data and test data

void RandomSamplingData(dMat& Ztrain, dMat& Ztest, dMat& Z);

void SamplingData(dMat& zTrain, dMat& Ztest, dMat& Z);

//---------------------------------------------------------------------
// cross-validation
void cvRandomSamplingData(dMat*& Ztrain, dMat*& Ztest, dMat& Z, char* filename);


#endif
