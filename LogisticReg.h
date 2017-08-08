#ifndef LOGISTIC_H_
#define LOGISTIC_H_

#include <NTL/ZZ.h>
#include <iostream>
#include <map>

#include "Cipher.h"
#include "Scheme.h"
#include "Database.h"



#include "LRtest.h"

using namespace std;
using namespace NTL;

class LR {
public:

    Scheme& scheme;
    LR(Scheme& scheme) : scheme(scheme) {};
    

    void EncryptData(Cipher*& zTrainCipher,  dMat zTrain, LRpar& LRparams, Scheme& scheme);
    
    void LogistRegHE(Cipher*& thetaCipher, Cipher*& zTrainCipher, LRpar& LRparams, Scheme& scheme, dMat zTrain);
    
    Cipher* getGrad_deg3(Cipher*& thetaCipher, Cipher*& zTrainCipher,  LRpar& LRparams);
    Cipher* getGrad_deg7(Cipher*& thetaCipher, Cipher*& zTrainCipher,  LRpar& LRparams);
    
    
};



#endif /* LOGISTIC_H_ */
