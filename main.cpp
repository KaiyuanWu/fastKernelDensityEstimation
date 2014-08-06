/* 
 * File:   main.cpp
 * Author: kaiwu
 *
 * Created on September 29, 2013, 9:38 PM
 */
#include <cstdlib>
#include <TProfile.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "iBalltree.h"
#include "TRandom.h"
#include "TFile.h"
#include "time.h"
#include "TPad.h"
#include "TH2F.h"
#include "TObject.h"
#include <string>
#include <TRandom1.h>
#include <TGraph.h>

#define PI 3.14159265353
using namespace std;

double hermitePolynomial(double x, int n){
    double h0=1.;
    double h1=2*x;
    double ret=0.;
    if(n==0)
        return h0;
    if(n==1)
        return h1;
    if(n>1){
        for(int i1=2;i1<=n;i1++){
            ret=2*x*h1-2*(i1-1)*h0;
            h0=h1;
            h1=ret;
        }
    }
    return ret;
}
double getGaussProbability(double* data,int* beta=NULL){
    double ret=1.;
    double H=1.;
    double mu=0.;
    double sigma=1.;
    for(int iDimension=0;iDimension<NDIMENSION;iDimension++){
        int b;
        if(!beta) b=0;
        else b=beta[iDimension];
        H=hermitePolynomial((data[iDimension]-mu)/(sigma*sqrt(2)),b)*pow(-1./(sigma*sqrt(2)),b);
        ret*=1./(sqrt(2*PI)*sigma)*exp(-0.5*pow((data[iDimension]-mu)/sigma,2.))*H;
    }
    return ret;
}

double testShuttle(int nPoints){
    TFile* fin=TFile::Open("exp/shuttle.scale.tr.root");
    int nDimension=NDIMENSION                  ;
    int  nTerms=8               ;
    int  nMinimumBoxSize=30     ;
    double epsilon=1.0e-2       ;
    double v[NDIMENSION];
    int label;
    TTree* tree=(TTree*) fin->Get("data");
    for(int iV=0;iV<NDIMENSION;iV++){
        tree->SetBranchAddress(Form("var%d",iV),&(v[iV]));
    }
    tree->SetBranchAddress("label",&label);
    int targetL=0;
    nPoints =0;
    int entries=tree->GetEntries();
    for(int iEntry=0;iEntry<entries;iEntry++){
        tree->GetEntry(iEntry);
        if(label==targetL)
            nPoints++;
    }
    double * data=new double[NDIMENSION*nPoints];
    double* weight=new double[nPoints]          ;
    int iPoint=0;
    for(int iEntry=0;iEntry<entries;iEntry++){
        tree->GetEntry(iEntry);
        if(label==targetL){
            for(int iDimension=0;iDimension<16;iDimension++)
                data[iPoint*NDIMENSION+iDimension]=v[iDimension];
            weight[iPoint]=1./nPoints;
            iPoint++;
        }
    }
    double h=.153            ;
    iBoxtree bT(nTerms, nMinimumBoxSize, epsilon) ;
    bT.setData(data,weight,h,nPoints)                    ;
    bT.buildTree()                                      ;
}
double testLetter(int nPoints){
    TFile* fin=TFile::Open("exp/letter.tr.root");
    int nDimension=NDIMENSION                  ;
    int* beta=new int[nDimension]       ;
    int  nTerms=8               ;
    int  nMinimumBoxSize=10     ;
    double epsilon=1.0e-2       ;
    for(int iDimension=0;iDimension<nDimension;iDimension++)
        beta[iDimension]=0; 
    double v[NDIMENSION];
    int label;
    TTree* tree=(TTree*) fin->Get("data");
    for(int iV=0;iV<NDIMENSION;iV++){
        tree->SetBranchAddress(Form("var%d",iV),&(v[iV]));
    }
    tree->SetBranchAddress("label",&label);
    int targetL=0;
    nPoints =0;
    int entries=tree->GetEntries();
    for(int iEntry=0;iEntry<entries;iEntry++){
        tree->GetEntry(iEntry);
        if(label==targetL)
            nPoints++;
    }
    double * data=new double[NDIMENSION*nPoints];
    double* weight=new double[nPoints]          ;
    int iPoint=0;
    for(int iEntry=0;iEntry<entries;iEntry++){
        tree->GetEntry(iEntry);
        if(label==targetL){
            for(int iDimension=0;iDimension<16;iDimension++)
                data[iPoint*NDIMENSION+iDimension]=v[iDimension];
            weight[iPoint]=1./nPoints;
            iPoint++;
        }
    }
    double h=.153            ;
    iBoxtree bT(nTerms, nMinimumBoxSize, epsilon) ;
    bT.setData(data,weight,h,nPoints)                    ;
    bT.buildTree()                                      ;
}
double testGauss(int nPoints){
    int nDimension=NDIMENSION                  ;
    int  nTerms=20                             ;
    int  nMinimumBoxSize= 10     ;
    double epsilon=1.0e-2       ;
    double* data=new double[nPoints*nDimension] ;
    double* weight=new double[nPoints]          ;
    double h=.15                               ;
    iBoxtree bT(nTerms, nMinimumBoxSize, epsilon) ;
    for(int iPoint=0;iPoint<nPoints;iPoint++){
        for(int iDimension=0;iDimension<nDimension;iDimension++)
            data[iPoint*nDimension+iDimension]=gRandom->Gaus(0,1.);
        weight[iPoint]=1./nPoints                       ;
    }
    bT.setData(data,weight,h,nPoints)                   ;
    bT.buildTree()                                      ;
    cout<<"here here!"<<endl;
    for(int iPoint=0;iPoint<10;iPoint++){
        cout<<bT.directEvaluate(data+iPoint*nDimension)<<" vs  "<<bT.fgtEvaluate(data+iPoint*nDimension)<<" "<<getGaussProbability(data+iPoint*nDimension)<<endl;
       //bT.fgtEvaluate(data+iPoint*nDimension);
    }
    return 0;
}
double testUniform(int nPoints){
    cout<<"nPoints= "<<nPoints<<endl;
    int nDimension=NDIMENSION                  ;
    int  nTerms=20              ;
    int  nMinimumBoxSize=log(nPoints)    ;
    double epsilon=1.0e-6       ;
    double* data=new double[nPoints*nDimension] ;
    double* weight=new double[nPoints]          ;
    double h=.17            ;
    iBoxtree bT(nTerms, nMinimumBoxSize, epsilon) ;
    for(int iPoint=0;iPoint<nPoints;iPoint++){
        for(int iDimension=0;iDimension<nDimension;iDimension++)
            data[iPoint*nDimension+iDimension]=gRandom->Rndm();
        weight[iPoint]=1./nPoints                       ;
    }
    bT.setData(data,weight,h,nPoints)                   ;
    bT.buildTree()                                      ;
    for(int iPoint=0;iPoint<nPoints;iPoint++){
         bT.fgtEvaluate(data+iPoint*nDimension);
//        //bT.directEvaluate(data+iPoint*nDimension);
        //cout<<bT.directEvaluate(data+iPoint*nDimension)<<" vs  "<<bT.fgtEvaluate(data+iPoint*nDimension)<<endl;
//        //cout<<bT.directEvaluate(data+iPoint*nDimension)<<endl;
    }
    return 0;
}
int main(int argc, char** argv) {
    gRandom->SetSeed()          ;
    int jobId=atoi(argv[1])     ;
    int nPoints=atoi(argv[2])   ;
    int nFunction=atoi(argv[3]) ;
    switch(nFunction){
        case 0:
            testShuttle(nPoints);
            break;
        case 1:
            testLetter(nPoints);
            break;
        case 2:
            testUniform(nPoints);
            break;
        case 3:
            testGauss(nPoints);
            break;
        default:
            cout<<"Error function type!"<<endl;
            exit(-1);
            break;
    }
    cout<<"nPoints= "<<nPoints<<endl;
    return 0;
}