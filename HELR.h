#ifndef LOGISTIC_H_
#define LOGISTIC_H_


#include <iostream>
#include <map>


#include <NTL/ZZ.h>
#include <NTL/BasicThreadPool.h>


#include "../src/TestScheme.h"
#include "../src/Cipher.h"
#include "../src/CZZ.h"
#include "../src/EvaluatorUtils.h"
#include "../src/Message.h"
#include "../src/Params.h"
#include "../src/PubKey.h"
#include "../src/Scheme.h"
#include "../src/SchemeAlgo.h"
#include "../src/SchemeAux.h"
#include "../src/SecKey.h"
#include "../src/StringUtils.h"



#include "Database.h"
#include "LRtest.h"

using namespace std;
using namespace NTL;


class LogReg {
public:

    Scheme& scheme;
    LRpar& LRparams;
    SecKey& secretKey;
    
    
    //! @ constructor
    LogReg(Scheme& scheme, SecKey& secretKey, LRpar& LRparams) : scheme(scheme),  secretKey(secretKey),  LRparams(LRparams) {}
    
    
    //---------------------------------------------------------------------------------------------------
    
    void EncryptData(Cipher*& zTrainCipher,  dMat zTrain);
    //void EncryptData_small(Cipher*& zTrainCipher,  dMat zTrain, LRpar& LRparams, Scheme& scheme);

    void HElogreg(Cipher*& thetaCipher, Cipher*& zTrainCipher,   dMat zTrain);
    
    Cipher* getgrad_deg3(Cipher*& thetaCipher, Cipher*& zTrainCipher);
    Cipher* getgrad_deg7(Cipher*& thetaCipher, Cipher*& zTrainCipher);
    
    static long getctlvl(Cipher& ctxt);
    void show_and_compare(dVec& theta, dMat zTrain, CZZ*& dtheta);

    
};



#endif /* LOGISTIC_H_ */
