/* 
 * File:   iBoxtree.h
 * Author: kaiwu
 *
 * Created on September 29, 2013, 9:39 PM
 */

#ifndef IBALLTREE_H
#define	IBALLTREE_H
#include <vector>
#include <iostream>
#include "math.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TBox.h"
#include <fstream>
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TNamed.h"
#include "TLegend.h"
#include "stdlib.h"
#include "TPad.h"
#include "TTree.h"
#include "TMarker.h"
#include "TPie.h"
#define NDIMENSION 1
#define SQRT2PI 2.50662827463100024
#define INVERSESQRT2PI 0.398942280401432703
#define SQRT2   1.41421356237309515
#define INVERSESQRT2 0.707106781186547462
using namespace std;

class iBoxtree {
public:
    iBoxtree(int nTerms, int nMinimumBoxSize, double epsilon);
    virtual ~iBoxtree();
    class iBoxtreeNode;
    //interaction list
    typedef struct INTERACTIVE_LIST {
        struct INTERACTIVE_LIST* _leftChild;
        struct INTERACTIVE_LIST* _rightChild;
        double _distance;
        iBoxtreeNode* _node;
    } INTERACTIVE_LIST_T;
    class iBoxtreeNode{
    public:
          iBoxtreeNode(iBoxtree* rootTree);
          ~iBoxtreeNode(){}          ;
          void calculateStatistics() ;
          void reset()               ;
          double fgtEvaluate(double* p)                 ;
          void   fgtEvaluate(double* p,int* beta,int totalNBetas,double* value)     ;
          void   calculateTargetCoefficients(double* p,int* nTerms) ;
          void   calculateTargetCoefficients(double* p) ;
          double calculateTargetError(double* p)        ;
          void   calculateTargetCoeffcients2(double*p)  ;
          int nPoints(){return _rightPoints-_leftPoints+1;}             ;
          void buildFGTCoefficient()                                    ;
          void buildFGTCoefficient2()                                   ;
          void printNode(const char* indent, bool last)                 ;
          void saveDiameter(vector<double> & diameter,vector<double> & diameterPCA,vector<double>& error);
          
          //statistics of this nodes
          double _mean[NDIMENSION]  ;
          double _sigma[NDIMENSION] ;
          double _min[NDIMENSION]   ;
          double _max[NDIMENSION]   ;
          int    _nTerms[NDIMENSION];
          int    _nTerms2[NDIMENSION];
          double _diameter[NDIMENSION]  ;
          double _diameterPCA[NDIMENSION]  ;
          double _weight        ;
          
          int _leftPoints       ;
          int _rightPoints      ;
          
          //variables of Fast Gauss Transform
          bool _hasFGTCoefficient       ;
          int _totalNTerms              ;
          int _totalNTerms2              ;
          
          double* _coefficients         ;
          double* _targetCoefficients   ;
          
          bool    _ableBeenApproxmiated ;
          
          double _error                 ;
                  
          iBoxtreeNode* _leftChild     ;
          iBoxtreeNode* _rightChild    ;
          iBoxtreeNode* _parent        ;
          iBoxtree*    _rootTree       ;
          
          double _pcaProjection[NDIMENSION*3];
          
          double _cut;
          int    _cutDimension;
          INTERACTIVE_LIST* _list      ;
    };
    
    
    
    //reset the interior attributes
    void buildTree()            ;   
    double hermitePolynomial(double x,int n);

    iBoxtreeNode* _root         ;
    void setData(double* data, double* weight,double bandwidth,int nPoint);
    double directEvaluate(double* p)    ;
    double fgtEvaluate(double* p)       ;
    double fgtEvaluate2(double* p)      ;
    //number of data points
    int _nPoints                ;
    double* _data               ;
    double* _weight             ;
    double _totalWeights        ;

    int _nMinimumBoxSize        ;
    int _nMaxTerms                 ;
    double _epsilon             ;
    //Inverse of bandwidth
    double _invBandwidth        ;
    double _cutOffDistance      ;
    double _scaleFactor         ;
    
private:
    //swap 2 data points
    void swap(int i,int j);
    //re-arrange the data points array
    void select(int dimension,int position, int low,int high);
    //build ball
    void buildBox(int low,int high,iBoxtreeNode* root);
    
    //parameters for the interpolation function
    double  _xH         ;
    double  _maxX       ;
    int     _minP       ;
    int     _maxP       ;
    int     _maxBeta    ;
    int     _nX         ;
    
    
};
#endif	/* IBALLTREE_H */