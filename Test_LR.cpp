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

    
    if(Argc != 2){
        cout << "-------------------------------------------------------------" << endl;
        cerr << "Enter the File and degree of approximation \t"  << "(e.g. $test edin.txt 3) \n ";
    }
    
    char* filename  =  Argv[1];
    
    
    
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
    long dim ;

    
    dMat mtheta3_list;
    dMat mtheta7_list;
    dMat mtheta_sig_list;
    
    ofstream fout;
    fout.open("beta.txt");

    
    for(int k = 0; k < 5; ++k){
        
        
        cout << "-------------------------------------------------------------" << endl;
        cout << "HELR_3 " << endl;
        
        int polydeg = 3;  // degree of approximation polynomial
        
        struct LRpar LRparams;
        ReadLRparams(LRparams, max_iter, zTrain[0], polydeg, logp);
        
        dVec mtheta3(LRparams.dim1, 0.0);
        for(int i= 0; i< LRparams.max_iter; i++){
            LR_poly(mtheta3, zTrain[k], LRparams);
        }
        
        for(int i= 0; i< LRparams.dim1; i++)
            cout << "[" << mtheta3[i] << "] " ;
        cout  << endl;
        
        getAUC(mtheta3, zTest[k]);
        mtheta3_list.push_back(mtheta3);
        
        
        cout << "-------------------------------------------------------------" << endl;
        cout << "HELR_7 " << endl;
        
        polydeg = 7;  // degree of approximation polynomial
        
        struct LRpar LRparams7;
        ReadLRparams(LRparams7, max_iter, zTrain[0], polydeg, logp);
        
        dVec mtheta7(LRparams7.dim1, 0.0);
        for(int i= 0; i< LRparams7.max_iter; i++){
            LR_poly(mtheta7, zTrain[k], LRparams7);
        }
        
        for(int i= 0; i< LRparams7.dim1; i++)
            cout << "[" << mtheta7[i] << "] " ;
        cout  << endl;
        
        getAUC(mtheta7, zTest[k]);
        mtheta7_list.push_back(mtheta7);
        
        dim = LRparams7.dim1;
        
        
        cout << "-------------------------------------------------------------" << endl;
        cout << "sigmoid LR " << endl;
        dVec mtheta_sig(LRparams.dim1, 0.0);
        
        for(int i= 0; i< LRparams.max_iter; i++){
            LR_sigmoid(mtheta_sig, zTrain[k], LRparams);
        }
        
        for(int i= 0; i< LRparams.dim1; i++)
            cout << "[" << mtheta_sig[i] << "] " ;
        cout  << endl;
        
        getAUC(mtheta_sig, zTest[k]);
        mtheta_sig_list.push_back(mtheta_sig);
        
        
        cout << "-------------------------------------------------------------" << endl;
        cout << "MSE (HELR/non-HELR): " << getMSE(mtheta3, mtheta_sig) << endl;
        cout << "MSE (HELR/non-HELR): " << getMSE(mtheta7, mtheta_sig) << endl;

    }
 
    
    //! write the beta results in the text file
    //fout << "-------------------------------------------------------------" << endl;
    //fout << "HELR_3" << endl;
    for(int i = 0; i < dim; ++i){
        for(int k = 0; k < 4; ++k){
            fout << mtheta3_list[k][i] << "," ;
        }
        fout << mtheta3_list[4][i] << ";" << endl;
    }
    
    
    //fout << "-------------------------------------------------------------" << endl;
    //fout << "HELR_7" << endl;
    for(int i = 0; i < dim; ++i){
        for(int k = 0; k < 4; ++k){
            fout << mtheta7_list[k][i] << "," ;
        }
        fout << mtheta7_list[4][i] << ";" << endl;
    }
    
    
    //fout << "-------------------------------------------------------------" << endl;
    //fout << "LR" << endl;
    for(int i = 0; i < dim; ++i){
        for(int k = 0; k < 4; ++k){
            fout << mtheta_sig_list[k][i] << "," ;
        }
        fout << mtheta_sig_list[4][i] << ";" << endl;
    }

    fout.close();
    
 
    
    delete[] zTest;
    delete[] zTrain;

    return 0;
}






