#include "LRtest.h"

ZZ ScaleUp(double value, const long& l){
    ZZ res;
    
    /*
    ZZ ztemp= to_ZZ("1");
    for(int i=0; i<l; i++)
        ztemp *= 2;
    
    res= to_ZZ(value*ztemp);
    */
    
    if(l < 31) {
        res= to_ZZ((1 << l) * value);
    }
    else {
        double dtemp= value * (1 << 30);
        double dtemp1= (1 << (l - 30));
        dtemp *= (dtemp1);
        
        res= to_ZZ(dtemp);
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
void LR_poly(dVec& theta, dMat zTrain, LRpar& LRparams){
    double* grad = new double[theta.size()];
    for(int i = 0; i < theta.size(); ++i) grad[i] = 0.0;
    
    // deg=3: (0.5- 1/2 (c1*ip+ c3*(ip)^3))z[k]
    for(int k = 0; k < zTrain.size(); ++k){
        double innerprod = inner_prod(zTrain[k], theta)/8;
        double dtemp= innerprod;
        
        double coeff = 0.0;
        
        for(int i= 1; i< 10; ++i){
            double power= pow(innerprod, i);
            coeff += LRparams.coeff[i] * power;
        }
        
        coeff= 0.5- coeff*(LRparams.polyscale)/8.0;
    
        for(int i = 0; i < theta.size(); ++i) grad[i] += coeff * zTrain[k][i];
    }
    
    for(int i = 0; i < theta.size(); ++i) theta[i] +=  grad[i]/zTrain.size();
    delete[] grad;
}




void getAUC(dVec theta, dMat zTest){
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
    
    cout << "Failure rate: (y = 1) " << n_fail_y1 << "/" << xtheta_y1.size() << " + (y = 0) " << n_fail_y0 << "/" ;
    cout << xtheta_y0.size() << " = " << 100.0 * (n_fail_y0 + n_fail_y1) / zTest.size() << " %." << endl;
    
    if(xtheta_y0.size() == 0 || xtheta_y1.size() ==0){ cout << "n_test_yi = 0 : cannot compute AUC" << endl; }
    else{
        double auc = 0.0;
        for(int i = 0; i < xtheta_y1.size(); ++i){ for(int j = 0; j < xtheta_y0.size(); ++j){
            if(xtheta_y0[j] <= xtheta_y1[i]) auc++; }}
        auc /= xtheta_y1.size() * xtheta_y0.size();
        cout << "AUC: " << auc << endl;
    }
}



// Update mtheta (in an unencrypted status) and compare the encrypted output "theta"
// with using the same sigmoid-approximation method

void show_and_compare(dVec& theta, dMat zTrain, CZZ*& dtheta, LRpar& LRparams, long logp){

    LR_poly(theta, zTrain, LRparams);
    
    double maxErrorBit = 0.0;
    double minRelativeBit= 20.0;
    
    for(int i=0; i< LRparams.dim1; i++){
        ZZ msg_scaleup= ScaleUp(theta[i], logp);
        
        cout << "m " << i << " : ["  <<  msg_scaleup << "], real theta: " << theta[i] << endl;   // unencrypted
        
        double theta;
        conv(theta, dtheta[i].r);
        for(int j=0; j< logp; j++){
            theta /= 2;
        }
        
        cout << "d " << i << " : ["  << dtheta[i].r << "], theta: " << theta << endl;      // encrypted
        cout << "e " << i << " : ["  <<  (msg_scaleup - dtheta[i].r)  << "],  ";    // error
        
        
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

