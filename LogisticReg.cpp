#include "LogisticReg.h"

#include <NTL/ZZ.h>
#include <cmath>
#include <map>
#include <math.h>  // pow
#include <sys/time.h>

#include "NumUtils.h"
#include "Params.h"
#include "PubKey.h"
#include "Scheme.h"
#include "SecKey.h"


#include "Ring2Utils.h"
#include "EvaluatorUtils.h"

#include <iostream>
#include <stdio.h>
#include <thread>

#include <vector>
#include <NTL/BasicThreadPool.h>
 
#include <sys/resource.h>   // check the memory usage
#include <sys/time.h>

#include "LRtest.h"




//! @ output ztrain = (2*y[k]-1) * (1, x[k])  -> encryption: zTrainCipher
void LR::EncryptData(Cipher*& zTrainCipher, dMat zTrain, LRpar& LRparams, Scheme& scheme){
    
    
    long slots = LRparams.n_training;  // same as n_training
    double* ctemp1= new double[LRparams.dim1];
    CZZ**  mvec= new CZZ*[LRparams.dim1];
    
    //1. zTrainCipher = (2*y[k]-1) * (1, x[k])  -> encryption: zTrainCipher
    for (long i = 0; i < LRparams.dim1; ++i){
        mvec[i]= new CZZ[slots];
        for(int k=0; k< slots; ++k){
            mvec[i][k] = EvaluatorUtils::evaluateVal(zTrain[k][i]/8, 0.0, scheme.params.logp);
        }
    }
    
    NTL_EXEC_RANGE(LRparams.dim1, first, last);
    for (long i = first; i < last; ++i){
        zTrainCipher[i] = scheme.encryptSlots(mvec[i], slots, scheme.params.logq);
    }
    NTL_EXEC_RANGE_END;

    
    delete[] ctemp1;
    delete[] mvec;

}




void LR::LogistRegHE(Cipher*& thetaCipher, Cipher*& zTrainCipher,  LRpar& LRparams, Scheme& scheme, dMat zTrain){
    struct rusage usage;
    auto start= chrono::steady_clock::now();
    
    // zSumCipher = 1/(n/polyscale) * sum (2y[k]-1) x[k] = 1/n * sum z[k]
    // At the first round, we just take zTrainCipher as theta (lvl= 2)
    Cipher* zSumCipher= new Cipher[LRparams.dim1];
    

    Cipher* ctemp= new Cipher[LRparams.dim1];
    long one= 1;
    
    NTL_EXEC_RANGE(LRparams.dim1, first, last);
    for (long i = first; i < last; ++i){
        zSumCipher[i]= zTrainCipher[i];
        
        for(long j=1; j< LRparams.logn; ++j){
            ctemp[i] = scheme.rotate2(zSumCipher[i], j);
            scheme.addAndEqual(zSumCipher[i], ctemp[i]);
        }
        ctemp[i] = scheme.rotate(zSumCipher[i], one);
        scheme.addAndEqual(zSumCipher[i], ctemp[i]);
        
        scheme.modSwitchAndEqual(zSumCipher[i], LRparams.logn - LRparams.log2polyscale);
        thetaCipher[i] = zSumCipher[i];
    }
    NTL_EXEC_RANGE_END;
    
    
    
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed= chrono::duration <double, milli> (diff).count()/1000.0;
    double totaltime = timeElapsed;
    
    int ret = getrusage(RUSAGE_SELF,&usage);
    cout << "-------------------------------------------------------------" << endl;
    cout << "1-iter : mod(theta)= " << thetaCipher[0].modulus << ", running time= "  << timeElapsed << "s, " ;
    cout<<  "Mem= " << usage.ru_maxrss/(1024)  << "MB" << endl;
    cout << "-------------------------------------------------------------" << endl;
    
    
    // compare and show between encrypted and unencrypted results
    dVec mtheta(LRparams.dim1, 0.0);
    CZZ* dtheta = new CZZ[LRparams.dim1];
    for(int i=0; i< LRparams.dim1; ++i){
        dtheta[i] = (scheme.decryptSlots(thetaCipher[i]))[0];
    }
    
    show_and_compare(mtheta, zTrain, dtheta, LRparams, scheme.params.logp);
    
    
    for(int j= 1; j< LRparams.max_iter; ++j){
        auto start= chrono::steady_clock::now();
        
        Cipher* gradCipher= new Cipher[LRparams.dim1];
        if(LRparams.polydeg==3){
            gradCipher= getGrad_deg3(thetaCipher, zTrainCipher, LRparams);
        }
        else{
            gradCipher= getGrad_deg7(thetaCipher, zTrainCipher, LRparams);
        }
    
        
        NTL_EXEC_RANGE(LRparams.dim1, first, last);
        for (long i = first; i < last; ++i){
            scheme.addAndEqual(thetaCipher[i], zSumCipher[i]);
            scheme.subAndEqual(thetaCipher[i], gradCipher[i]);
        }
        NTL_EXEC_RANGE_END;
        
        auto end = std::chrono::steady_clock::now();
        auto diff = end - start;
        double timeElapsed= chrono::duration <double, milli> (diff).count()/1000.0;
        totaltime += timeElapsed;
        
        ret = getrusage(RUSAGE_SELF,&usage);
        
        delete[] gradCipher;
        
        cout << "-------------------------------------------------------------" << endl;
        cout << j+1 << "-iter : mod(theta)= " << thetaCipher[0].modulus << ", running time= "  << timeElapsed << "s, " ;
        cout<<  "Mem= " << usage.ru_maxrss/(1024)  << "MB" << endl;
        cout << "-------------------------------------------------------------" << endl;
        
        
        for(int i=0; i< LRparams.dim1; ++i){
            dtheta[i] = (scheme.decryptSlots(thetaCipher[i]))[0];
        }
        
        show_and_compare(mtheta, zTrain, dtheta, LRparams, scheme.params.logp);
    }
    
    cout << "Total Evaluation Time = " << totaltime << " s" << endl;

    delete[] dtheta;
}



 
Cipher* LR::getGrad_deg7(Cipher*& thetaCipher, Cipher*& zTrainCipher, LRpar& LRparams){
    
    //! compute ztheta =  (p* z[k]/8) (p* theta)/p  : mod(theta)+ logp
    long modtheta= thetaCipher[0].modulus;
    Cipher* ctemp= new Cipher[LRparams.dim1];
    
    NTL_EXEC_RANGE(LRparams.dim1, first, last);
    for (long i = first; i < last; i++){
        ctemp[i] = zTrainCipher[i];
        scheme.multAndEqual(ctemp[i], thetaCipher[i]);
    }
    NTL_EXEC_RANGE_END;
    
    Cipher ztheta = ctemp[0];
    for (long i = 1; i < LRparams.dim1; i++){
        scheme.addAndEqual(ztheta, ctemp[i]);
    }
     
    scheme.modSwitchAndEqual(ztheta, scheme.params.logp);

    //-----------------------------------------------------------------------------------------
    // ctemp  =  alpha * (ztheta) * (a7*ztheta^6 + a5*ztheta^4 +a3*ztheta^2 + a1) * z[i] for 0<=i<=dim
    //        =  ((alpha*a7*z[i]) * (ztheta)) * (ztheta^6 + a5/a7*ztheta^4 +a3/a7*ztheta^2 + a1/a7)
    //-----------------------------------------------------------------------------------------
   

    // zSquare = p * (theat*z[k]/8)^2
    Cipher zSquare= scheme.square(ztheta);
    scheme.modSwitchAndEqual(zSquare, scheme.params.logp);

    // zQuartic = p * (theat*z[k]/8)^4
    Cipher zQuartic = scheme.square(zSquare);
    scheme.modSwitchAndEqual(zQuartic, scheme.params.logp);

    // zQuartic = p * (a7*ztheta^4 + a5*ztheta^2 + a3)
    Cipher ctemp1= scheme.multByConst(zSquare, LRparams.evalcoeff[5]);
    scheme.modSwitchAndEqual(ctemp1, scheme.params.logp);
    scheme.addConstAndEqual(ctemp1, LRparams.evalcoeff[3]);
    
    scheme.multByConst(zQuartic, LRparams.evalcoeff[7]);
    scheme.addAndEqual(zQuartic, ctemp1);
   

    Cipher* res= new Cipher[LRparams.dim1];
    
    NTL_EXEC_RANGE(LRparams.dim1, first, last);
    for (long i = first; i < last; ++i){
        res[i] = zTrainCipher[i];   
        scheme.multAndEqual(res[i], ztheta);      
        scheme.modSwitchAndEqual(res[i], scheme.params.logp);  
    
        ctemp[i]= scheme.multByConst(res[i], LRparams.evalcoeff[1]);
        scheme.modSwitchAndEqual(ctemp[i], scheme.params.logp);
    
        scheme.multAndEqual(res[i], zSquare);     
        scheme.modSwitchAndEqual(res[i], scheme.params.logp);  
        

        scheme.multAndEqual(res[i], zQuartic);     
        scheme.modSwitchAndEqual(res[i], scheme.params.logp);  
        scheme.multByConst(res[i], LRparams.evalcoeff[7]);

        scheme.addAndEqual(res[i], ctemp[i]);
    }
    NTL_EXEC_RANGE_END;
    
     
    //! AllSum
    long one= 1;
    NTL_EXEC_RANGE(LRparams.dim1, first, last);
    for (long i = first; i < last; ++i){
        for(long j=1; j< LRparams.logn; ++j){
            ctemp[i] = scheme.rotate2(res[i], j);
            scheme.addAndEqual(res[i], ctemp[i]);
        }
        ctemp[i] = scheme.rotate(res[i], one);
        scheme.addAndEqual(res[i], ctemp[i]);

        scheme.modSwitchAndEqual(res[i], LRparams.logn - LRparams.log2polyscale);
    }
    NTL_EXEC_RANGE_END;

    delete[] ctemp;
    return res;
}



//-----------------------------------------------------------------------------------------
// res[i]  =  1/n * sum_k (ztheta) * (a3*ztheta^3 + a1*ztheta) * z[i]/8 for 0<=i<=dim
//         =  1/n * sum_k ((z[i]/8) * (ztheta)) * (ztheta^2 + a1/a3) * a3)
//-----------------------------------------------------------------------------------------

Cipher* LR::getGrad_deg3(Cipher*& thetaCipher, Cipher*& zTrainCipher, LRpar& LRparams){
    
    //! compute ztheta =  (p* z[k]/8) (p* theta)/p  : mod(theta)+ logp
    long modtheta= thetaCipher[0].modulus;
    Cipher* ctemp= new Cipher[LRparams.dim1];
    
    NTL_EXEC_RANGE(LRparams.dim1, first, last);
    for (long i = first; i < last; i++){
        ctemp[i] = zTrainCipher[i];
        scheme.multAndEqual(ctemp[i], thetaCipher[i]);
    }
    NTL_EXEC_RANGE_END;
    
    Cipher ztheta = ctemp[0];
    for (long i = 1; i < LRparams.dim1; i++){
        scheme.addAndEqual(ztheta, ctemp[i]);
    }
    
    scheme.modSwitchAndEqual(ztheta, scheme.params.logp);
    
    
    //! zSquare= p * (theat*z[k]/8)^2 + p (a1/a3)
    Cipher zSquare= scheme.square(ztheta);
    scheme.modSwitchAndEqual(zSquare, scheme.params.logp);
    scheme.addConstAndEqual(zSquare, LRparams.evalcoeff[1]);
    
    
    //! res[i]= ((z[i]/8) * (ztheta)) * (ztheta^2 + a1/a3) * a3)
    Cipher* res= new Cipher[LRparams.dim1];
    
    NTL_EXEC_RANGE(LRparams.dim1, first, last);
    for (long i = first; i < last; ++i){
        res[i]= scheme.multByConst(zTrainCipher[i], LRparams.evalcoeff[3]);
        scheme.modSwitchAndEqual(res[i], scheme.params.logp);
        
        scheme.multAndEqual(res[i], ztheta);
        scheme.modSwitchAndEqual(res[i], scheme.params.logp);
        
        scheme.multAndEqual(res[i], zSquare);
        scheme.modSwitchAndEqual(res[i], scheme.params.logp);
    }
    NTL_EXEC_RANGE_END;
    
    
    //! (sum res[i]) / (n/polyscale)
    long one= 1;
    
    NTL_EXEC_RANGE(LRparams.dim1, first, last);
    for (long i = first; i < last; ++i){
        for(long j=1; j< LRparams.logn; ++j){
            ctemp[i] = scheme.rotate2(res[i], j);
            scheme.addAndEqual(res[i], ctemp[i]);
        }
        ctemp[i] = scheme.rotate(res[i], one);
        scheme.addAndEqual(res[i], ctemp[i]);
        
        scheme.modSwitchAndEqual(res[i], LRparams.logn -  LRparams.log2polyscale);
    }
    NTL_EXEC_RANGE_END;
    
    delete[] ctemp;
    return res;
}

 

