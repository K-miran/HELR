#ifndef LR_UNENC
#define LR_UNENC

#include <NTL/ZZ.h>
#include "Database.h"
#include "CZZ.h"

using namespace std;
using namespace NTL;

ZZ ScaleUp(double value, const long& l);


//---------------------------------------------------------------------
//! @ Plaintext logistic regression
void LR_sigmoid(dVec& theta, dMat zTrain, LRpar& LRparams);

//! @ Linear approximation
void LR_linear(dVec& theta, dMat zTrain, LRpar& LRparams);

void LR_cubic(dVec& theta, dMat zTrain, LRpar& LRparams);

void LR_poly(dVec& theta, dMat zTrain, LRpar& LRparams);

void show_and_compare(dVec& theta, dMat zTrain, CZZ*& dtheta, LRpar& LRparams, long logp);


void getAUC(dVec theta, dMat zTest);

#endif
