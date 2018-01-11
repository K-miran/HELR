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


//---------------------------------------------------------------------
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

//---------------------------------------------------------------------

void SamplingData(dMat& Ztrain, dMat& Ztest, dMat& Z){
    int n_train = 1 << ((int) log2(Z.size() * 0.91));
    for(int k = 0; k < n_train; ++k) Ztrain.push_back(Z[k]);
    for(int k = 0; k < Z.size() - n_train; ++k) Ztest.push_back(Z[k]);
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


//---------------------------------------------------------------------

void cvRandomSamplingData(dMat*& Ztrain, dMat*& Ztest, dMat& Z, char* filename){
    
    int n_train;
    int n_test[5];
    
    string str(filename);
     
    if(str == "data/edin.txt"){
        n_train  = 1024;
        n_test[0]= 250;
        n_test[1]= 250;
        n_test[2]= 251;
        n_test[3]= 251;
        n_test[4]= 251;
    }
    else if(str == "data/lbw.txt"){
        n_train =  256;
        n_test[0]= 37;
        n_test[1]= 38;
        n_test[2]= 38;
        n_test[3]= 38;
        n_test[4]= 38;
    }
    else if(str == "data/nhanes3.txt"){
        n_train = 16384;
        n_test[0]= 3129;
        n_test[1]= 3130;
        n_test[2]= 3130;
        n_test[3]= 3130;
        n_test[4]= 3130;
    }
    else if(str == "data/pcs.txt"){
        n_train =  512;
        n_test[0]= 75;
        n_test[1]= 76;
        n_test[2]= 76;
        n_test[3]= 76;
        n_test[4]= 76;
    }
    else if(str == "data/uis.txt"){
        n_train =  512;
        n_test[0]= 115;
        n_test[1]= 115;
        n_test[2]= 115;
        n_test[3]= 115;
        n_test[4]= 115;
    }
    
   
    
    
    // Z0: zero vector of dimension "d"
    vector<double> Z0;
    for(int l = 0; l < Z[0].size(); ++l){
        Z0.push_back(0);
    }
    
    vector<bool> randBit(Z.size(), false);
    srand(time(NULL));
    
    //--------------------------------------------------
    //! Generate the testing dataset
    int k;
    ofstream fout;
    fout.open("test_data.txt");

    
    for(int l = 0; l < 4; ++l){
        k = 0;
    
        while(k < n_test[l]){
            int i = rand() % Z.size();
            if(randBit[i] == false){
                randBit[i] = true;
                Ztest[l].push_back(Z[i]);
                for(int l = 0; l < Z[i].size(); ++l){
                    fout << Z[i][l]<< ",";
                }
                fout << endl;
                
                k++;
            }
        }
        //fout << "----------------------------------------------" << endl;
    }
    
    
    
    for(int i = 0; i < Z.size(); ++i){
        if(randBit[i] == false){
            Ztest[4].push_back(Z[i]);
            for(int l = 0; l < Z[i].size(); ++l){
               fout << Z[i][l]<< ",";
            }
            fout << endl;
            
        }
    }
    fout.close();
    
    
    //--------------------------------------------------
    //! Generate the learning dataset : Ztrain[0] ... Ztrain[4]
    //! Ztrain[0] = (Ztest[1], Ztest[2], Ztest[3], Ztest[4])
    //! Ztrain[1] = (Ztest[2], Ztest[3], Ztest[4], Ztest[0])
    
    for(int m = 0; m < 5; ++m){
        for(int l = m+1; l < 5; ++l){
            for(int i = 0; i < n_test[l]; ++i){
                Ztrain[m].push_back(Ztest[l][i]);
                //cout << Ztest[l][i][0] << "\t" ;
            }
            //cout << endl;
        }

        for(int l = 0; l < m; ++l){
            for(int i = 0; i < n_test[l]; ++i){
                Ztrain[m].push_back(Ztest[l][i]);
                //cout << Ztest[l][i][0] << "\t" ;
            }
            //cout << endl;
        }
        
        int l1 = n_train  + n_test[m] - Z.size();
        
        for(int l = 0; l< l1; ++l){
            Ztrain[m].push_back(Z0);
        }
        
        cout << "#(learning samples[" << m << "]): " << Z.size() -  n_test[m] << ", " ;
        cout << "Ztrain.size: " << Ztrain[m].size() << endl;
    }
}








