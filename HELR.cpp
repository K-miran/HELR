
#include <cmath>
#include <map>
#include <math.h>  // pow
#include <sys/time.h>
#include <sys/resource.h>   // check the memory usage

#include <iostream>
#include <stdio.h>
#include <thread>
#include <vector>

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


#include "LRtest.h"
#include "HELR.h"

#define debug 0


long LogReg::getctlvl(Cipher& ctxt){
    long nBits = NumBits(ctxt.mod)-1;
    return nBits;
}



//---------------------------------------------------------------------------------------------------

//! @ output ztrain = (2*y[k]-1) * (1, x[k])  -> encryption: zTrainCipher
void LogReg::EncryptData(Cipher*& zTrainCipher, dMat zTrain){
    
    //long slots = LRparams.n_training;
    
    CZZ**  mvec= new CZZ*[LRparams.dim1];
    
    // zTrainCipher = (2*y[k]-1) * (1, x[k])
    for (long i = 0; i < LRparams.dim1; ++i){
        mvec[i]= new CZZ[LRparams.nslots];
        for(int k=0; k< LRparams.nslots; ++k){
            if(zTrain[k][i]!=0){
                mvec[i][k].r = scaleup(zTrain[k][i], LRparams.logp-3);  // logp/8
            }
            else{
                mvec[i][k].r = to_ZZ("0");
            }
            mvec[i][k].i = to_ZZ("0");
        }
    }
    
    NTL_EXEC_RANGE(LRparams.dim1, first, last);
    for (long i = first; i < last; ++i){
        zTrainCipher[i] = scheme.encrypt(mvec[i], LRparams.nslots);
    }
    NTL_EXEC_RANGE_END;

   
    delete[] mvec;

}


void LogReg::show_and_compare(dVec& theta, dMat zTrain, CZZ*& dtheta){
    
    LR_poly(theta, zTrain, LRparams);
    
    double maxErrorBit = 0.0;
    double minRelativeBit= 20.0;
    
    for(int i=0; i< LRparams.dim1; i++){
        ZZ msg_scaleup= scaleup(theta[i], LRparams.logp);
        
        cout << "m " << i << " : ["  <<  msg_scaleup << "], theta: " << theta[i] << endl;   // unencrypted
        
        double theta;
        conv(theta, dtheta[i].r);
        theta = scaledown(theta, LRparams.logp);
        
        cout << "d " << i << " : ["  << dtheta[i].r << "], (HE)theta: " << theta << endl;      // encrypted
        cout << "e " << i << " : ["  << (msg_scaleup - dtheta[i].r)  << "],  ";    // error
        
        
        // Msg-bit
        ZZ msgbnd= dtheta[i].r;
        double MsgBit = 0.0;
        if(msgbnd!=0)
            MsgBit= (log(abs(msgbnd))/log(2)) ;
        cout << "Msg: " << MsgBit << ", ";
        
        // Error-bit
        ZZ error = (msg_scaleup - dtheta[i].r);
        double ErrorBit = 0.0;
        if(error!=0)
            ErrorBit= (log(abs(error))/log(2)) ;
        cout << "Error: " << ErrorBit << endl;
        
        double RelativeBit= MsgBit- ErrorBit;
        
        if(maxErrorBit < ErrorBit) maxErrorBit= ErrorBit;
        if(minRelativeBit > RelativeBit)  minRelativeBit= RelativeBit;
        cout << "-------------------------------------------------------------" << endl;
        
    }
    
    cout << "MAX error bit : " << maxErrorBit << ", Min relative-error bit: " << minRelativeBit<< endl;
    cout << "-------------------------------------------------------------" << endl;
 
    
}



void LogReg::HElogreg(Cipher*& thetaCipher, Cipher*& zTrainCipher,  dMat zTrain){
    struct rusage usage;
    auto start= chrono::steady_clock::now();
    
    // zSumCipher = 1/(n/polyscale) * sum (2y[k]-1) x[k] = 1/n * sum z[k]
    // At the first round, we just take zTrainCipher as theta (lvl= 2)
    Cipher* zSumCipher= new Cipher[LRparams.dim1];
    
    //! rot and sum
    Cipher* ctemp= new Cipher[LRparams.dim1];
    
    
    NTL_EXEC_RANGE(LRparams.dim1, first, last);
    for (long i = first; i < last; ++i){
        zSumCipher[i]= zTrainCipher[i];
        
        for(long j= 0; j< LRparams.logn; ++j){
            int l = (1<<j);
            ctemp[i] = scheme.leftRotate(zSumCipher[i], l);
            scheme.addAndEqual(zSumCipher[i], ctemp[i]);
        }
    
        
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
    cout << "1-iter : mod(theta)= " << getctlvl(thetaCipher[0]) << ", running time= "  << timeElapsed << "s, " ;
    cout<<  "Mem= " << usage.ru_maxrss/(1024*1024)  << "GB" << endl;
    cout << "-------------------------------------------------------------" << endl;
    
    

    dVec mtheta(LRparams.dim1, 0.0);
    CZZ* dtheta = new CZZ[LRparams.dim1];
    for(int i = 0; i< LRparams.dim1; ++i){
        dtheta[i] = (scheme.decrypt(secretKey, thetaCipher[i]))[0];
    }
    
    show_and_compare(mtheta, zTrain, dtheta);
 
    
    //--------------------------------------------------------------
    double zlvl =  getctlvl(zSumCipher[0]);
    
    for(int j= 1; j< LRparams.max_iter; ++j){
        auto start= chrono::steady_clock::now();
        
        Cipher* gradCipher= new Cipher[LRparams.dim1];
        
        switch(LRparams.polydeg){
            case 3:
                gradCipher= getgrad_deg3(thetaCipher, zTrainCipher);
                break;
            case 7:
                gradCipher= getgrad_deg7(thetaCipher, zTrainCipher);
                break;
                
        }
        
        double tlvl = getctlvl(thetaCipher[0]);
        double glvl = getctlvl(gradCipher[0]);
        
        NTL_EXEC_RANGE(LRparams.dim1, first, last);
        for (long i = first; i < last; ++i){
            Cipher ztemp = scheme.modEmbed(zSumCipher[i], zlvl- tlvl);
            scheme.addAndEqual(thetaCipher[i], ztemp);
            
            scheme.modEmbedAndEqual(thetaCipher[i], tlvl- glvl);
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
        cout << j+1 << "-iter : mod(theta)= " << getctlvl(thetaCipher[0]) << ", running time= "  << timeElapsed << "s, " ;
        cout<<  "Mem= " << usage.ru_maxrss/(1024*1024)  << "GB" << endl;
        cout << "-------------------------------------------------------------" << endl;
        
        
        for(int i=0; i< LRparams.dim1; ++i){
            dtheta[i] = (scheme.decrypt(secretKey, thetaCipher[i]))[0];
        }
        
        show_and_compare(mtheta, zTrain, dtheta);
    }
    
    
    cout << "Total Evaluation Time = " << totaltime << " s" << endl;

    delete[] zSumCipher;
    delete[] ctemp;
    delete[] dtheta;
   
    
}



 
Cipher* LogReg::getgrad_deg7(Cipher*& thetaCipher, Cipher*& zTrainCipher){
    
    //! compute ztheta =  (p* z[k]/8) (p* theta)/p  : mod(theta)+ logp
    double tlvl =  getctlvl(thetaCipher[0]);
    double zlvl =  getctlvl(zTrainCipher[0]);
    
    Cipher* ctemp= new Cipher[LRparams.dim1];
    
    NTL_EXEC_RANGE(LRparams.dim1, first, last);
    for (long i = first; i < last; i++){
        ctemp[i] = scheme.modEmbed(zTrainCipher[i], zlvl - tlvl);
        scheme.multAndEqual(ctemp[i], thetaCipher[i]);
    }
    NTL_EXEC_RANGE_END;
    
    Cipher ztheta = ctemp[0];
    for (long i = 1; i < LRparams.dim1; i++){
        scheme.addAndEqual(ztheta, ctemp[i]);
    }
     
    scheme.modSwitchAndEqual(ztheta, LRparams.logp);  // ztheta.lvl : tlvl - logp

    
    
    //-----------------------------------------------------------------------------------------
    // ctemp  =  alpha * (ztheta) * (a7*ztheta^6 + a5*ztheta^4 +a3*ztheta^2 + a1) * z[i] for 0<=i<=dim
    //        =  ((alpha*a7*z[i]) * (ztheta)) * (ztheta^6 + a5/a7*ztheta^4 +a3/a7*ztheta^2 + a1/a7)
    //-----------------------------------------------------------------------------------------
   

    // zSquare = p * (theat*z[k]/8)^2  with  tlvl - 2*logp
    Cipher zSquare= scheme.square(ztheta);
    scheme.modSwitchAndEqual(zSquare, LRparams.logp);

    // zQuartic = p * (theat*z[k]/8)^4 with  tlvl - 3*logp
    Cipher zQuartic = scheme.square(zSquare);
    scheme.modSwitchAndEqual(zQuartic, LRparams.logp);

    // zQuartic = p * (a7*ztheta^4 + a5*ztheta^2 + a3)  with  tlvl - 3*logp
    Cipher ctemp1= scheme.multByConst(zSquare, LRparams.evalcoeff[5]);
    scheme.modSwitchAndEqual(ctemp1, LRparams.logp);     // lvl(theta)+3
    scheme.addConstAndEqual(ctemp1, LRparams.evalcoeff[3]);
    
    scheme.multByConst(zQuartic, LRparams.evalcoeff[7]);
    scheme.addAndEqual(zQuartic, ctemp1);
   

    Cipher* res= new Cipher[LRparams.dim1];
    
    NTL_EXEC_RANGE(LRparams.dim1, first, last);
    for (long i = first; i < last; ++i){
        res[i] = scheme.modEmbed(zTrainCipher[i], zlvl - tlvl + LRparams.logp);
        scheme.multAndEqual(res[i], ztheta);      
        scheme.modSwitchAndEqual(res[i], LRparams.logp);       // res: tlvl - logp
    
        ctemp[i]= scheme.multByConst(res[i], LRparams.evalcoeff[1]);
        scheme.modSwitchAndEqual(ctemp[i], LRparams.logp);     // res: tlvl - 2*logp
    
        scheme.multAndEqual(res[i], zSquare);     
        scheme.modSwitchAndEqual(res[i], LRparams.logp);      // res: tlvl - 3 *logp
        

        scheme.multAndEqual(res[i], zQuartic);     
        scheme.modSwitchAndEqual(res[i], LRparams.logp);
        scheme.multByConst(res[i], LRparams.evalcoeff[7]);    // res: tlvl - 4 *logp

        scheme.addAndEqual(res[i], ctemp[i]);
    }
    NTL_EXEC_RANGE_END;
    
     
    //! rot and sum
   
    NTL_EXEC_RANGE(LRparams.dim1, first, last);
    for (long i = first; i < last; ++i){
        for(long j= 0; j< LRparams.logn; ++j){
            long l = (1<<j);
            ctemp[i] = scheme.leftRotate(res[i], l);
            scheme.addAndEqual(res[i], ctemp[i]);
        }
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

Cipher* LogReg::getgrad_deg3(Cipher*& thetaCipher, Cipher*& zTrainCipher){
    
    double tlvl =  getctlvl(thetaCipher[0]);
    double zlvl =  getctlvl(zTrainCipher[0]);
    
    Cipher* ctemp= new Cipher[LRparams.dim1];
    
    NTL_EXEC_RANGE(LRparams.dim1, first, last);
    for (long i = first; i < last; i++){
        ctemp[i] = scheme.modEmbed(zTrainCipher[i], zlvl - tlvl);
        scheme.multAndEqual(ctemp[i], thetaCipher[i]);
    }
    NTL_EXEC_RANGE_END;
    
    Cipher ztheta = ctemp[0];
    for (long i = 1; i < LRparams.dim1; i++){
        scheme.addAndEqual(ztheta, ctemp[i]);
    }
    
    scheme.modSwitchAndEqual(ztheta, LRparams.logp);  // ztheta.lvl : tlvl - logp

 
    //---------------------------------------------------
    //! zSquare= p * (theat*z[k]/8)^2 + p (a1/a3) with  tlvl - 2*logp
    Cipher zSquare= scheme.square(ztheta);
    scheme.modSwitchAndEqual(zSquare, LRparams.logp);
    scheme.addConstAndEqual(zSquare, LRparams.evalcoeff[1]);
    
 

    //! res[i]= ((z[i]/8) * a3 ) * (ztheta)) * (ztheta^2 + a1/a3))
    Cipher* res= new Cipher[LRparams.dim1];
    
    NTL_EXEC_RANGE(LRparams.dim1, first, last);
    for (long i = first; i < last; ++i){
        res[i]= scheme.multByConst(zTrainCipher[i], LRparams.evalcoeff[3]);
        scheme.modSwitchAndEqual(res[i], LRparams.logp);  // ((z[i]/8) * a3 ) with zlvl - logp
        
        scheme.modEmbedAndEqual(res[i], zlvl - tlvl);
        scheme.multAndEqual(res[i], ztheta);
        scheme.modSwitchAndEqual(res[i], LRparams.logp);  // ((z[i]/8) * a3 ) * (ztheta)) with zlvl - 2*logp
        
        scheme.multAndEqual(res[i], zSquare);
        scheme.modSwitchAndEqual(res[i], LRparams.logp);  // with zlvl - 3*logp
    }
    NTL_EXEC_RANGE_END;
  
    
    
    //! (sum res[i]) / (n/polyscale)

    NTL_EXEC_RANGE(LRparams.dim1, first, last);
    for (long i = first; i < last; ++i){
        for(long j= 0; j< LRparams.logn; ++j){
            long l = (1<<j);
            ctemp[i] = scheme.leftRotate(res[i], l);
            scheme.addAndEqual(res[i], ctemp[i]);
        }
       
        scheme.modSwitchAndEqual(res[i], LRparams.logn - LRparams.log2polyscale);
    }
    NTL_EXEC_RANGE_END;
    
    delete[] ctemp;
    return res;
}


