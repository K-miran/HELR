#include <cmath>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <sys/time.h>
#include <chrono>
#include <sys/resource.h>   // check the memory usage
#include <stdio.h>
#include <thread>
#include <fstream>
#include <sstream>

#include <NTL/RR.h>
#include <NTL/xdouble.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/BasicThreadPool.h>


#include "../src/CZZ.h"
#include "../src/Params.h"
#include "../src/PubKey.h"
#include "../src/Scheme.h"
#include "../src/SchemeAlgo.h"
#include "../src/SecKey.h"
#include "../src/TestScheme.h"
#include "../src/TimeUtils.h"
#include "../src/Ring2Utils.h"
#include "../src/StringUtils.h"
#include "../src/EvaluatorUtils.h"

#include "Database.h"
#include "LRtest.h"
#include "HELR.h"


using namespace NTL;
using namespace std;



int main(int Argc, char** Argv) {

    
    if(Argc != 3){
        cout << "-------------------------------------------------------------" << endl;
        cerr << "Enter the File and degree of approximation \t"  << "(e.g. $test edin.txt 3) \n ";
    }
    
    char* filename  =  Argv[1];
    int polydeg = atoi(Argv[2]);  // degree of approximation polynomial
    
    
    dMat  zData;
    dMat* zTest = new dMat[5];
    dMat* zTrain = new dMat[5];
    
    
    int nLine= readData(zData, filename);
    
    cout << "Sample the learning and test data ..." << endl;
    cvRandomSamplingData(zTrain, zTest, zData, filename);

    
    //----------------------------------------------------------------
    // Parameters for Logistic regression
    //----------------------------------------------------------------

    long logN= 11;
    long logp= 28;
    long logl= 10;
    long logq, cBit1, cBit2;
    int max_iter;
    
    struct LRpar LRparams;
    ReadLRparams(LRparams, max_iter, zTrain[0], polydeg, logp);
    
    SetNumThreads(LRparams.dim1);
    //SetNumThreads(4);
    
    switch(polydeg){
        case 3:
            cBit1=  (LRparams.logn - LRparams.log2polyscale);          // 1st iteration
            cBit2 =  (3*logp+ LRparams. logn - LRparams.log2polyscale);  // 2nd~ iteration
            logq = cBit1 + (LRparams.max_iter-1)*(cBit2)+ logp + logl;                  // max-bitlength we need
            break;
            
        case 7:
            cBit1=  (LRparams.logn - LRparams.log2polyscale);          // 1st iteration
            cBit2=  (4*logp+ LRparams.logn - LRparams.log2polyscale);
            logq= cBit1 + (LRparams.max_iter-1)*(cBit2)+ logp + logl;
            break;
    }
    
 
    
    
    cout << "Data dimension with dummy vectors: " << LRparams.dim1 << ", Number of lines: " << nLine << endl;
    
    cout << "-------------------------------------------------------------" << endl;
    cout << "Key Generation ... (logN,logp,logq, nslots)= (" ;
    cout << logN << "," << logp << "," << logq << "," << LRparams.nslots << ")" <<endl;
    
    auto start= chrono::steady_clock::now();
    
    Params params(logN, logq);
    SecKey secretKey(params);
    PubKey publicKey(params, secretKey);
    SchemeAux schemeaux(logN);
    Scheme scheme(params, publicKey, schemeaux);
    SchemeAlgo algo(scheme);
    
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    cout << "KeyGen time= " << chrono::duration <double, milli> (diff).count()/1000.0 << " s" << endl;

    

    
    LogReg LR(scheme, secretKey, LRparams);
    dMat HEtheta_list;
    
    ofstream fout;
    fout.open("beta.txt");

    
    for(int k = 0; k < 1; ++k){
        
        cout << "-------------------------------------------------------------" << endl;
        cout << k << "th Data Encryption ... " << endl;
        
        struct rusage usage;
        
        start= chrono::steady_clock::now();
        
        Cipher* zTrainCipher = new Cipher[LRparams.dim1];
        
        LR.EncryptData(zTrainCipher, zTrain[k]);
        
        end = std::chrono::steady_clock::now();
        diff = end - start;
        cout << "Enc time= "  << chrono::duration <double, milli> (diff).count()/1000.0 << "(s), " ;
        
        int ret = getrusage(RUSAGE_SELF,&usage);
        cout<< "Mem: " << usage.ru_maxrss/(1024)  << "(MB)" << endl;
        
    
        cout << "-------------------------------------------------------------" << endl;
        cout << "HE Logistic Regression ... "  << endl;
        
        Cipher* thetaCipher= new Cipher[LRparams.dim1];
        
        start= chrono::steady_clock::now();
        
        LR.HElogreg(thetaCipher, zTrainCipher, zTrain[k]);
        
        
        end = std::chrono::steady_clock::now();
        diff = end - start;
        cout << "Eval time= "  << chrono::duration <double, milli> (diff).count()/1000.0 << "(s), " ;
        
        ret = getrusage(RUSAGE_SELF,&usage);
        cout<< "Mem: " << usage.ru_maxrss/(1024)  << "(MB)" << endl;
        
        
        cout << "-------------------------------------------------------------" << endl;
        cout << "Decryption ... "  << endl;
        
        dVec HEtheta(LRparams.dim1, 0.0);

        CZZ* dtheta = new CZZ[LRparams.dim1];
    
        for(int i=0; i< LRparams.dim1; ++i){
            dtheta[i] = (scheme.decrypt(secretKey, thetaCipher[i]))[0];
            
            conv(HEtheta[i], dtheta[i].r);
            HEtheta[i] = scaledown(HEtheta[i], LRparams.logp);
            cout << "[" << HEtheta[i] << "] " ;
        }
        cout << ": enc " << endl;
        
        getAUC(HEtheta, zTest[k]);
        HEtheta_list.push_back(HEtheta);
        
        
        cout << "-------------------------------------------------------------" << endl;
        cout << "Compare with unenc LR " << endl;
        dVec mtheta(LRparams.dim1, 0.0);
        
        for(int i= 0; i< LRparams.max_iter; i++){
            LR_poly(mtheta, zTrain[k], LRparams);
        }
        
        for(int i= 0; i< LRparams.dim1; i++)
            cout << "[" << mtheta[i] << "] " ;
        cout << ": unenc " << endl;
        
        getAUC(mtheta, zTest[k]);
        
        cout << "MSE (HELR/non-HELR): " << getMSE(HEtheta, mtheta) << endl;
        
        
        cout << "-------------------------------------------------------------" << endl;
        cout << "Compare with sigmoid LR " << endl;
        dVec mtheta_sig(LRparams.dim1, 0.0);
        
        for(int i= 0; i< LRparams.max_iter; i++){
            LR_sigmoid(mtheta_sig, zTrain[k], LRparams);
        }
        
        for(int i= 0; i< LRparams.dim1; i++)
            cout << "[" << mtheta_sig[i] << "] " ;
        cout << ": unenc " << endl;
        
        getAUC(mtheta_sig, zTest[k]);
        
        cout << "MSE (HELR/non-HELR): " << getMSE(HEtheta, mtheta_sig) << endl;

    }
 
    
    //! write the beta results in the text file
    fout << "-------------------------------------------------------------" << endl;
    fout << "[" << endl;
    for(int i = 0; i < LRparams.dim1; ++i){
        for(int k = 0; k < 4; ++k){
            fout << HEtheta_list[k][i] << "," ;
        }
        fout << HEtheta_list[4][i] << ";" << endl;
    }
    fout << "];" << endl;
    fout.close();
    
 
    
    delete[] zTest;
    delete[] zTrain;

    return 0;
}






