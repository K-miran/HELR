#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "math.h"

#include "time.h"
#include "Database.h"

using namespace std;



int readData(dMat& Z, char* filename){
    string str(filename);
    //    ifstream infile("data/edin.txt");
    ifstream infile(str);
    string line;
    char split_char = '\t';
    int dim1;
    
    if(getline(infile, line)){
        istringstream split(line);
        vector<string> tokens;
        for (string each; getline(split, each, split_char); tokens.push_back(each));
        dim1 = tokens.size();
        
        vector<double> Z0;
        Z0.push_back(stod(tokens[dim1 - 1]) * 2 - 1);
        if(Z0[0] != 1 && Z0[0] != -1){ cout << "Error: Y[0] is non-binary" << endl; }
        for(int i = 0; i < dim1 - 1; ++i) Z0.push_back(stod(tokens[i]) * Z0[0]);
        Z.push_back(Z0);
    }
    else{ cout << "Error: file reading error" << endl; }
    
    while(getline(infile, line)){
        istringstream split(line);
        vector<string> tokens;
        for (string each; getline(split, each, split_char); tokens.push_back(each));
        if(tokens.size() != dim1){ cout << "Error: database dimension" << endl; }
        vector<double> Zi;
        Zi.push_back(stod(tokens[dim1 - 1]) * 2 - 1);
        if(Zi[0] != 1 && Zi[0] != -1){ cout << "Error: Y[i] is non-binary" << endl; }
        for(int i = 0; i < dim1 - 1; ++i) Zi.push_back(stod(tokens[i]) * Zi[0]);
        Z.push_back(Zi);
    }
    
    int nLine = Z.size();
    
    double XjMax, XjMean;
    
    for(int j = 1; j < dim1; j++){
        XjMax = 0.0;
        for(int i = 0; i < nLine; ++i){
            if(XjMax < fabs(Z[i][j])) XjMax = fabs(Z[i][j]);
        }
        if(XjMax > 1){
            XjMean = 0.0;
            for(int i = 0; i < nLine; ++i) XjMean += fabs(Z[i][j]);
            XjMean /= nLine;
            for(int i = 0; i < nLine; ++i) Z[i][j] /= XjMean;
        }
    }
    return nLine;
}

void RandomSamplingData(dMat& Ztrain, dMat& Ztest, dMat& Z){
    int n_train = 1 << ((int) log2(Z.size() * 0.91));
    vector<bool> randBit(Z.size(), false);
    srand(time(NULL));
    
    int k = 0;
    while(k < n_train){
        int i = rand() % Z.size();
        if(randBit[i] == false){
            randBit[i] = true;
            Ztrain.push_back(Z[i]);
            k++;
        }
    }
    
    for(int i = 0; i < Z.size(); ++i){
        if(randBit[i] == false){ Ztest.push_back(Z[i]); }
    }
}

void SamplingData(dMat& Ztrain, dMat& Ztest, dMat& Z){
    int n_train = 1 << ((int) log2(Z.size() * 0.91));
    for(int k = 0; k < n_train; ++k) Ztrain.push_back(Z[k]);
    for(int k = 0; k < Z.size() - n_train; ++k) Ztest.push_back(Z[k]);
}



void ReadLRparams(LRpar& LRparams, int max_iter, dMat zTrain, double* coeff, ZZ* evalcoeff, int polyscale, int polydeg, long Slots){
    LRparams.max_iter= max_iter;
    
    LRparams.dim1 = zTrain[0].size();
    LRparams.n_training = zTrain.size();
    LRparams.logn = log2(LRparams.n_training);
    
    for(int i=0; i<10; i++){
        LRparams.coeff[i] = coeff[i];
        LRparams.evalcoeff[i] = evalcoeff[i];
    }
    
    LRparams.polyscale= polyscale;
    LRparams.log2polyscale= log2(polyscale);
    
    LRparams.polydeg= polydeg;
    
    LRparams.Slots= Slots;

}
