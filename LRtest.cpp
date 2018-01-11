#include "LRtest.h"




void ReadLRparams(LRpar& LRparams, int max_iter, dMat zTrain,  int polydeg,  long logp){
    LRparams.max_iter= max_iter;
    
    LRparams.dim1 = zTrain[0].size();
    //LRparams.n_training = zTrain.size();
    LRparams.n_training = zTrain.size();   // Power of two
    LRparams.logn = log2(LRparams.n_training);
    
    
    //----------------------------------------------
    // polynomial coefficients
    for(int i=0; i<10; i++){
         LRparams.coeff[i] = 0.0;
    }
    
    // "Scaled" polynomial coefficients
    LRparams.polyscale= 4;
    
    for(int i=0; i<10; i++){
        LRparams.evalcoeff[i]= to_ZZ("0");
    }

    
    switch(polydeg){
        case 3:
            LRparams.max_iter = 25;
            LRparams.coeff[1] =  2.401926612;
            LRparams.coeff[3] = (-1.631249824);
            
            LRparams.evalcoeff[1]=  scaleup(LRparams.coeff[1]/LRparams.coeff[3], logp);
            LRparams.evalcoeff[3]=  scaleup(LRparams.coeff[3], logp);
            break;
            
        case 7:
            LRparams.max_iter =   20;
            LRparams.coeff[1] =   3.46992;
            LRparams.coeff[3] = - 8.38814;
            LRparams.coeff[5] =   10.86804;
            LRparams.coeff[7] = - 5.0;
            
            LRparams.evalcoeff[1]=  scaleup(LRparams.coeff[1], logp);
            LRparams.evalcoeff[3]=  scaleup(LRparams.coeff[3], logp);
            LRparams.evalcoeff[5]=  scaleup(LRparams.coeff[5], logp);
            LRparams.evalcoeff[7]=  to_ZZ(LRparams.coeff[7]);
            break;
    }

    
    
    LRparams.log2polyscale= log2(LRparams.polyscale);
    LRparams.polydeg= polydeg;
    
    LRparams.logp = logp;
    
    LRparams.nslots= zTrain.size();  // POT close to "N*0.8"
    
}


//---------------------------------------------------------------------


ZZ scaleup(double value, const long& l){
    ZZ precision = power2_ZZ(l);
    double dtemp;
    conv(dtemp, precision);
    dtemp *= value;
    ZZ res = to_ZZ(dtemp);
    
    return res;
}



double scaledown(double value, const long& l){
    double res = value;
    for(int j=0; j< l; j++){
        res /= 2;
    }
    return res;
    
}



double inner_prod(dVec u, dVec v, int start = 0){
    double res = 0.0;
    for(int i = start; i < u.size(); ++i) res += u[i] * v[i];
    return res;
}

//---------------------------------------------------------------------

//! One iteration with the original method
void LR_sigmoid(dVec& theta, dMat zTrain, LRpar& LRparams){
    double* grad = new double[theta.size()];
    for(int i = 0; i < theta.size(); ++i) grad[i] = 0.0;
    for(int k = 0; k < zTrain.size(); ++k){
        double coeff = 1.0 / (1.0 + exp(inner_prod(zTrain[k], theta)));
        for(int i = 0; i < theta.size(); ++i) grad[i] += coeff * zTrain[k][i];
    }
    for(int i = 0; i < theta.size(); ++i) theta[i] += grad[i]/zTrain.size();;
    delete[] grad;
}

 

// update for one iteration
// deg3: 0.5 -
// deg7: 0.5 - (3.4 (x/8)+ ... )
// 4/n sum_i  (1- \sum coeff[j] * (ztrain*beta/8)^j ) z[i]/8
// 1/n sum_i  (0.5 - 0.5 \sum coeff[j] * (ztrain*beta/8)^j ) z[i]
void LR_poly(dVec& theta, dMat zTrain, LRpar& LRparams){
    double* grad = new double[theta.size()];
    for(int i = 0; i < theta.size(); ++i) grad[i] = 0.0;
    
    
    // deg=3: (0.5- 1/2 (c1*ip+ c3*(ip)^3))z[k]
    for(int k = 0; k < zTrain.size(); ++k){
        double innerprod = inner_prod(zTrain[k], theta)/8.0;
        double dtemp= innerprod;
        
        double coeff = 0.0;
        
        for(int i= 1; i< 10; ++i){
            double power= pow(innerprod, i);
            coeff += LRparams.coeff[i] * power;
        }
        
        coeff= 1.0 - coeff;
        
        for(int i = 0; i < theta.size(); ++i) grad[i] += coeff * (zTrain[k][i]/8.0);
    }
    
    for(int i = 0; i < theta.size(); ++i) theta[i] +=  grad[i]/(zTrain.size()/4.0);
    delete[] grad;
}




double getAUC(dVec theta, dMat zTest){
    int n_fail_y1 = 0;
    int n_fail_y0 = 0;
    
    dVec xtheta_y1;
    dVec xtheta_y0;
    
    for(int i = 0; i < zTest.size(); ++i){
        if(zTest[i][0] == 1.0){
            if(inner_prod(zTest[i], theta) < 0) n_fail_y1++;
            xtheta_y1.push_back(zTest[i][0] * inner_prod(zTest[i], theta, 1));
        }
        else{
            if(inner_prod(zTest[i], theta) < 0) n_fail_y0++;
            xtheta_y0.push_back(zTest[i][0] * inner_prod(zTest[i], theta, 1));
        }
    }
    
    double correctness= 100.0 - (100.0 * (n_fail_y0 + n_fail_y1) / zTest.size());
    cout << "Failure rate: (y = 1) " << n_fail_y1 << "/" << xtheta_y1.size() << " + (y = 0) " << n_fail_y0 << "/" ;
    cout << xtheta_y0.size() << " = " <<  (100.0 * (n_fail_y0 + n_fail_y1) / zTest.size()) << " %." << endl;
    cout << "Correctness: " << correctness  << " %." << endl;
    
    
    if(xtheta_y0.size() == 0 || xtheta_y1.size() ==0){
        cout << "n_test_yi = 0 : cannot compute AUC" << endl;
        return 0.0;
    }
    else{
        double auc = 0.0;
        for(int i = 0; i < xtheta_y1.size(); ++i){
            for(int j = 0; j < xtheta_y0.size(); ++j){
                if(xtheta_y0[j] <= xtheta_y1[i]) auc++;
            }
        }
        auc /= xtheta_y1.size() * xtheta_y0.size();
        return auc;
        cout << "AUC: " << auc << endl;
    }
}

double  getMSE(dVec theta1, dVec theta2){
    double res= 0.0;
    
    for(int i=0; i< theta1.size(); ++i){
        double dtemp = pow(theta1[i]-theta2[i], 2.0);
        res += dtemp;
    }
    res/= (theta1.size());
    return res;
    
}


double  getNMSE(dVec theta1, dVec theta2){
    double res= 0.0;
    
    for(int i=0; i< theta1.size(); ++i){
        double dtemp = pow(theta1[i], 2.0);
        res += dtemp;
    }
    res/= (theta1.size());
    
    double mse= getMSE(theta1, theta2);
    res= (mse/res);
    
    return res;
    
}



// Update mtheta (in an unencrypted status) and compare the encrypted output "theta"
// with using the same sigmoid-approximation method

//void show_and_compare(dVec& theta, dMat zTrain, CZZ*& dtheta, LRpar& LRparams, long logp){

