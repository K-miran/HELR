#include "CZZ.h"
#include "Params.h"
#include "PubKey.h"
#include "Scheme.h"
#include "SecKey.h"


#include "Ring2Utils.h"
#include "EvaluatorUtils.h"
#include "Database.h"
#include "LRtest.h"
#include "LogisticReg.h"

#include "NTL/ZZ.h"
#include "NTL/ZZX.h"
#include <NTL/BasicThreadPool.h>


#include <iostream>
#include <stdio.h>
#include <thread>
#include <vector>

#include <sys/resource.h>   // check the memory usage
#include <sys/time.h>


#include <cmath>
#include <chrono>

using namespace NTL;
using namespace std;




int main(int Argc, char** Argv) {

    if(Argc != 3){
        cout << "-------------------------------------------------------------" << endl;
        cerr << "Enter the File and degree of approximation \t"  << "(e.g. $test edin.txt 3) \n ";
    }
 
    cout << "-------------------------------------------------------------" << endl;
    cout << "Read the data ... " << endl;
    
    char* filename  =  Argv[1];
    int polydeg = atoi(Argv[2]);
    
    dMat Z, zTest, zTrain;

    int nLine= readData(Z, filename);
    
    cout << "Sample the learning and test data ..." << endl;
    
    //SamplingData(zTrain, zTest, Z);
    RandomSamplingData(zTrain, zTest, Z);
    
    
    //----------------------------------------------------------------
    // Parameters for Logistic regression
    //----------------------------------------------------------------
    int max_iter;          // number of iterations
    
    if(polydeg==3){
        max_iter= 25;
    }
    else{
        max_iter= 20;
    }
    
    long logN= 11;          // ring dimension
    long logp= 28;          // scale factor
    long logl= 10;
    long logq0= logp+ logl; // the smallest modulus of ciphertext

    
    // coefficients of polynomial for evaluation of GD (scaled by "2")
    double* coeff= new double[10];
    ZZ* evalcoeff= new ZZ[10];
    
    for(int i=0; i<10; i++){
        coeff[i]= 0.0;
        evalcoeff[i]= to_ZZ("0");
    }


    if(polydeg==3){
        coeff[1]=  2.401926612;
        coeff[3]= (-1.631249824);
        
        evalcoeff[1]=  ScaleUp(coeff[1]/coeff[3], logp);
        evalcoeff[3]=  ScaleUp(coeff[3], logp);
    }
    else{
        coeff[1]=   3.46992;
        coeff[3]= - 8.38814;
        coeff[5]=   10.86804;
        coeff[7]= - 5.0;
        
        evalcoeff[1]=  ScaleUp(coeff[1], logp);
        evalcoeff[3]=  ScaleUp(coeff[3], logp);
        evalcoeff[5]=  ScaleUp(coeff[5], logp);
        evalcoeff[7]=  to_ZZ("-5");
    }
    
    int polyscale= 4;       // 2 (scale) * polyscale = "8"
 
    
   
    
    long slots= zTrain.size();
    
    struct LRpar LRparams;
    ReadLRparams(LRparams, max_iter, zTrain, coeff, evalcoeff, polyscale, polydeg, slots);
    
    
    
    long consumed_mod1, consumed_mod2, logq;
    
    consumed_mod1=  (LRparams.logn - LRparams.log2polyscale);              // 1st iteration
    
    if(polydeg==3){
        consumed_mod2=  (3*logp+ LRparams.logn - LRparams.log2polyscale);  // 2nd~ iteration
        logq= consumed_mod1 + (max_iter-1)*(consumed_mod2)+ logq0;         // max-bitlength we need
        cout << "Deg 3-approximation with " ;
    }
    else{
        consumed_mod1=  (LRparams.logn - LRparams.log2polyscale);
        consumed_mod2=  (4*logp+ LRparams.logn - LRparams.log2polyscale);
        logq= consumed_mod1 + (max_iter-1)*(consumed_mod2)+ logq0;
        cout << "Deg 7-approximation with " ;
    }
    
    
    SetNumThreads(LRparams.dim1);
    
    cout << " Data dimension (X,Y): " << LRparams.dim1 << ", Number of lines: " << nLine ;
    cout << ", Number of learning samples: " << LRparams.n_training <<endl;
 

    cout << "-------------------------------------------------------------" << endl;
    cout << "Key Generation ... (logN,logp,logqL)= (" ;
    cout << logN << "," << logp << "," << logq << ")" <<endl;

    auto start= chrono::steady_clock::now();
    
    Params params(logN, logp, logq);
    SecKey secretKey(params);
    PubKey publicKey(params, secretKey);
    
    SchemeAux schemeaux(logp, logN + 2);
    Scheme scheme(params, secretKey, publicKey, schemeaux);

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    cout << chrono::duration <double, milli> (diff).count()/1000.0 << " s" << endl;
    
    LR SGD(scheme);

 
    cout << "-------------------------------------------------------------" << endl;
    cout << "Data Encryption ... " ;
    

    struct rusage usage;
    
    start= chrono::steady_clock::now();
    
    Cipher* zTrainCipher = new Cipher[LRparams.dim1];

    SGD.EncryptData(zTrainCipher, zTrain, LRparams, scheme);

    end = std::chrono::steady_clock::now();
    diff = end - start;
    cout << chrono::duration <double, milli> (diff).count()/1000.0 << "(s), " ;

    int ret = getrusage(RUSAGE_SELF,&usage);
    cout<< "Mem: " << usage.ru_maxrss/(1024)  << "(MB)" << endl;
    
 
    cout << "-------------------------------------------------------------" << endl;
    cout << "Evaluation ... " << endl;

    Cipher* thetaCipher= new Cipher[LRparams.dim1];
 
    SGD.LogistRegHE(thetaCipher, zTrainCipher,  LRparams, scheme, zTrain);
 
 

    cout << "-------------------------------------------------------------" << endl;
    cout << "-------------------------------------------------------------" << endl;
    
    cout << "Compare with unenc/enc approximate method ... " << endl;
    dVec theta_poly(LRparams.dim1, 0.0);
    
    for(int i= 0; i< LRparams.max_iter; i++)
        LR_poly(theta_poly, zTrain, LRparams);
    
    for(int i= 0; i< LRparams.dim1; i++)
        cout << "[" << theta_poly[i] << "] " ;
    cout << ": unenc " << endl;
    
    
 
    dVec theta_new(LRparams.dim1);
    for(int i= 0; i< LRparams.dim1; i++){
        CZZ* dvec = scheme.decryptSlots(thetaCipher[i]);
        conv(theta_new[i], dvec[0].r);
        for(int j=0; j< scheme.params.logp; j++) theta_new[i] /= 2;
    }
    
    for(int i=0; i< LRparams.dim1; i++) cout << "[" << theta_new[i] << "] " ;
    cout << ": enc " << endl;

    
    getAUC(theta_new, zTest);
 
  

    return 0;
}






