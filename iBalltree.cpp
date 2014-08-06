/* 
 * File:   iBoxtree.cpp
 * Author: kaiwu
 * 
 * Created on September 29, 2013, 9:39 PM
 */

#include "iBalltree.h"
#include "initParameters.h"
#include <iomanip>
#include <TGraph.h>
#include <bitset>
#include <list>
//initialization
//input parameters:derivative array bets
//                 number of terms
//                 minimum box size
//this array is served as a  temporary array
double targetCoefficients[1024]  ;
iBoxtree::iBoxtree(int nMaxTerms, int nMinimumBoxSize, double epsilon) {
    //parameters for the error estimation
    _nPoints = 0        ;
    _xH      = 0.01     ;
    _minP    = 2        ;
    _maxP    = 29       ;
    _maxBeta = 20       ;
    _maxX    = 4.       ;
    _nX      = 400      ;
    //check number of truncated terms
    _nMaxTerms = nMaxTerms;
    if (_nMaxTerms < _minP) {
        cout << "nTerms " << _nMaxTerms << " should greater than " << _minP << endl;
        _nMaxTerms = _minP;
    } else if (_nMaxTerms > _maxP) {
        cout << "nTerms " << _nMaxTerms << " should less than " << _maxP << endl;
        _nMaxTerms = _maxP;
    }
    //check the minimum box size
    _nMinimumBoxSize = nMinimumBoxSize;
    _epsilon = epsilon;
    //initialization no tree structure, no data, no weight
    _root   = NULL      ;
    _data   = NULL      ;
    _weight = NULL      ;
    _epsilon=epsilon    ;
}
//De-construction function 
iBoxtree::~iBoxtree() {
    if (_root != NULL) {
        _root->reset();
        delete _root;
    }
    if (_data != NULL)
        delete[] _data;
    if (_weight != NULL)
        delete[] _weight;
}
//set data
void iBoxtree::setData(double* data, double* weight,double bandwidth,int nPoint) {
    _nPoints = nPoint;
    //If the data once has been initialized, we need to clear the memory
    if (_data) {
        cout << "data has been set! reinitialize the data!" << endl;
        delete[] _data;
    }
    if (_weight)
        delete[] _weight;
    _data = new double[nPoint * NDIMENSION];
    memcpy(_data, data, _nPoints*NDIMENSION*sizeof (double));
    _weight = new double[nPoint];
    memcpy(_weight, weight, _nPoints * sizeof (double));
    if(bandwidth<=0){
        cout<<"Error bandwidth "<<bandwidth<<", Please check!"<<endl;
        exit(-1);
    }
    _totalWeights=0.;
    for(int iPoint=0;iPoint<nPoint;iPoint++)
        _totalWeights+=_weight[iPoint];
    _invBandwidth=1./bandwidth;
    _scaleFactor =1.          ;
    for(int iDimension=0;iDimension<NDIMENSION;iDimension++)
        _scaleFactor*=(INVERSESQRT2PI*_invBandwidth);
    //cut off distance
    // -h**2*log(epsilon (pi)**d/2 h**d)
    _epsilon/=fabs(_scaleFactor);
    _cutOffDistance=-2*bandwidth*bandwidth*log(_epsilon);
}

//initialize the tree Node
iBoxtree::iBoxtreeNode::iBoxtreeNode(iBoxtree* rootTree) {
    _leftPoints = _rightPoints = 0      ;
    _leftChild  = _rightChild  = NULL   ;
    _weight     = 0.                    ;
    _rootTree   = rootTree              ;
    _totalNTerms= 0                     ;
    _hasFGTCoefficient = false          ;
    _targetCoefficients=NULL            ;
    _coefficients=NULL                  ;
    _error=1.0e300                      ;
}
void iBoxtree::iBoxtreeNode::calculateStatistics() {
    _weight = 0.;
    for (int iDimension = 0; iDimension < NDIMENSION; iDimension++) {
        _mean[iDimension]=0.            ;
        _sigma[iDimension]=0.           ;
        _min[iDimension]=1.0e300        ;
        _max[iDimension]=-1.0e300       ;
    }
    for (int iPoint = _leftPoints; iPoint <= _rightPoints; iPoint++) {
        for (int iDimension = 0; iDimension < NDIMENSION; iDimension++) {
            _mean[iDimension]  += _rootTree->_data[iPoint * NDIMENSION + iDimension] * _rootTree->_weight[iPoint];
            _sigma[iDimension] += _rootTree->_data[iPoint * NDIMENSION + iDimension] * _rootTree->_data[iPoint * NDIMENSION + iDimension] * _rootTree->_weight[iPoint];
            if (_min[iDimension] > _rootTree->_data[iPoint * NDIMENSION + iDimension]) _min[iDimension] = _rootTree->_data[iPoint * NDIMENSION + iDimension];
            if (_max[iDimension] < _rootTree->_data[iPoint * NDIMENSION + iDimension]) _max[iDimension] = _rootTree->_data[iPoint * NDIMENSION + iDimension];
        }
        _weight += _rootTree->_weight[iPoint];
    }
    for (int iDimension = 0; iDimension < NDIMENSION; iDimension++) {
        _mean[iDimension]  /= _weight;
        _sigma[iDimension] /= _weight;
        _sigma[iDimension] -= _mean[iDimension] * _mean[iDimension]     ;
        _sigma[iDimension]  = sqrt(fabs(_sigma[iDimension]))            ;
        double r=_max[iDimension]-_mean[iDimension]                     ;
        double s=_mean[iDimension]-_min[iDimension]                     ;
        _diameter[iDimension]=r<s? s:r;
    }
}
void iBoxtree::iBoxtreeNode::reset() {
    if (_leftChild != NULL) {
        _leftChild->reset();
        delete _leftChild;
    }
    if (_rightChild != NULL) {
        _rightChild->reset();
        delete _rightChild;
    }
    if(_targetCoefficients) delete[] _targetCoefficients;
    if(_coefficients)       delete[] _coefficients;
}

void iBoxtree::buildTree() {
    //if previously has built a tree, we first delete this tree
    if (_root != NULL){
        _root->reset();
        delete _root;
    }
    //initialize the building...
    _root = new iBoxtreeNode(this);
    _root ->_parent=NULL;
    if (_totalWeights != 1.) {
        for (int i1 = 0; i1 < _nPoints; i1++) {
            _weight[i1] /= _totalWeights;
        }
    }
    buildBox(0, _nPoints - 1, _root)      ;
    
//    _root->_list=(INTERACTIVE_LIST* )malloc(sizeof(INTERACTIVE_LIST));
//    _root->_list->_node=_root            ;
//    _root->_list->_leftChild=NULL        ;
//    _root->_list->_rightChild=NULL       ;
//    calculateList(_root->_leftChild )    ;
//    calculateList(_root->_rightChild)    ;
//    _root->printNo    de("",0)              ;  
    _root->buildFGTCoefficient()        ;
//    ///////////////////////////////////////////////////
//    _root->printNode("",0)              ;    
    vector<double> diameters    ;
    vector<double> diametersPCA ;
    vector<double> error;
    _root->saveDiameter(diameters,diametersPCA,error);
    TFile* fout=TFile::Open("saveDiameter.root","update");
    TH1F* herror=new TH1F("herror","herror",1000,-10,10);
    for (int iDimension = 0; iDimension < NDIMENSION; iDimension++) {
        TH1F* h1 = new TH1F(Form("h diameter%d", iDimension), Form("h dimension :%d", iDimension), 1000, 0, 10);
        TH1F* h2 = new TH1F(Form("h pca diameter%d", iDimension), Form("h pca dimension :%d", iDimension), 1000, 0, 10);
        for (int is = 0; is < diameters.size()/NDIMENSION; is++){
            h1->Fill(diameters[is*NDIMENSION+iDimension]);
            h2->Fill(diametersPCA[is*NDIMENSION+iDimension]);
        }
        h1->Write("", TObject::kOverwrite);
        h2->Write("", TObject::kOverwrite);
    }
    for (int is = 0; is < diameters.size()/NDIMENSION; is++){
        if(error[is]==0.)
            herror->Fill(0.);
        else
            herror->Fill(log(error[is]));
    }
    herror->Write("",TObject::kOverwrite);
    fout->Close();
//    ////////////////////////////////////////////////////
}
void iBoxtree::iBoxtreeNode::saveDiameter(vector<double> & diameter,vector<double> & diameterPCA,vector<double>& error){
    if (!_leftChild&&!_rightChild){
        for (int iDimension = 0; iDimension < NDIMENSION; iDimension++){
            diameter.push_back(_diameter[iDimension] * _rootTree->_invBandwidth);
            diameterPCA.push_back(_diameterPCA[iDimension] * _rootTree->_invBandwidth);
        }
        error.push_back(_error);
    }
    if(_leftChild)
        _leftChild->saveDiameter(diameter,diameterPCA,error);
    if(_rightChild)
        _rightChild->saveDiameter(diameter,diameterPCA,error);
}
void iBoxtree::iBoxtreeNode::printNode(const char* indent, bool last){
//    if(!_leftChild&&!_rightChild){
//        printf("Node (%d,%d)\n",_leftPoints,_rightPoints);
//        _rootTree->printList("",_list,false);
//    }
//    else{
//        if(_leftChild)
//            _leftChild->printNode("",false);
//        if(_rightChild)
//            _rightChild->printNode("",false);
//    }
//    return;
    bool printDetail=true;
    if (_leftChild&&_rightChild) {
        if(last){
            printf("%s|-(%d,%d) ",indent,_leftPoints,_rightPoints);
            if (printDetail) {
                for (int iDimension = 0; iDimension < NDIMENSION; iDimension++)
                    printf(" %1.3f ", _diameter[iDimension]);
                for (int iDimension = 0; iDimension < NDIMENSION; iDimension++)
                    printf(" [%1.3f,%1.3f, %1.3f] ", _min[iDimension], _mean[iDimension], _max[iDimension]);
            }
            printf(" Error= %f \n",_error);
            _leftChild->printNode(Form("%s|  ",indent), true);
            _rightChild->printNode(Form("%s|  ",indent), false);
        }
        else{
            printf("%s+-(%d,%d) ",indent,_leftPoints,_rightPoints);
            if (printDetail) {
                for (int iDimension = 0; iDimension < NDIMENSION; iDimension++)
                    printf(" %1.3f ", _diameter[iDimension]);
                for (int iDimension = 0; iDimension < NDIMENSION; iDimension++)
                    printf(" [%1.3f,%1.3f, %1.3f] ", _min[iDimension], _mean[iDimension], _max[iDimension]);
            }
            printf(" Error= %f \n",_error);
            _leftChild ->printNode(Form("%s  ",indent), true) ;
            _rightChild->printNode(Form("%s  ",indent), false);
        }
    } else{
        if(!last){
            printf("%s+-(%d,%d) ",indent,_leftPoints,_rightPoints);
            if (printDetail) {
                for (int iDimension = 0; iDimension < NDIMENSION; iDimension++)
                    printf(" %1.3f ", _diameter[iDimension]);
                for (int iDimension = 0; iDimension < NDIMENSION; iDimension++)
                    printf(" [%1.3f,%1.3f, %1.3f] ", _min[iDimension], _mean[iDimension], _max[iDimension]);
            }
            printf(" Error= %f \n",_error);
        }
        else{
            printf("%s|-(%d,%d) ",indent,_leftPoints,_rightPoints);
            if (printDetail) {
                for (int iDimension = 0; iDimension < NDIMENSION; iDimension++)
                    printf(" %1.3f ", _diameter[iDimension]);
                for (int iDimension = 0; iDimension < NDIMENSION; iDimension++)
                    printf(" [%1.3f,%1.3f, %1.3f] ", _min[iDimension], _mean[iDimension], _max[iDimension]);
            }
            printf(" Error= %f \n",_error);
        }
//        for(int iP=_leftPoints;iP<=_rightPoints;iP++){
//            printf("(%f,%f) ",_rootTree->_data[iP*NDIMENSION],_rootTree->_data[iP*NDIMENSION+1]);
//        }
//        printf("\n");
    }
}
double iBoxtree::iBoxtreeNode::fgtEvaluate(double* p) {
    double leftValue = 0., rightValue = 0.;
    double ret = 0.;
    //check whether this node can be directly eliminated
    double minimumDistance=0.;
    for(int iDimension=0;iDimension<NDIMENSION;iDimension++){
        double r=0.;
        if(p[iDimension]<_min[iDimension])
            r=_min[iDimension]-p[iDimension];
        else if(p[iDimension]>_max[iDimension])
            r=p[iDimension]-_max[iDimension];
        else
            r=0.;
        minimumDistance+=r*r;
    }
    if(minimumDistance>_rootTree->_cutOffDistance)
        return 0.;
    //Terminal Node
    if (_leftChild == NULL && _rightChild == NULL) {
        if (_hasFGTCoefficient) {
            double targetError = calculateTargetError(p);
            if (targetError * _error < _rootTree->_epsilon) {
                calculateTargetCoefficients(p);
                ret = 0.;
                for (int iTerm = 0; iTerm < _totalNTerms; iTerm++){
                    ret += _coefficients[iTerm] * targetCoefficients[iTerm];
                }
                ret *= targetError*targetError;
                return ret;
            }
        }
        //we have to directly calculate the summation
        for (int iPoint = _leftPoints; iPoint <= _rightPoints; iPoint++) {
            double temp2 = 0.;
            for (int iDimension = 0; iDimension < NDIMENSION; iDimension++) {
                double r=(p[iDimension] - _rootTree->_data[iPoint*NDIMENSION+iDimension]) * _rootTree->_invBandwidth;
                temp2  += r*r;
            }
            temp2 *= (-0.5);
            ret += exp(temp2) * _rootTree->_weight[iPoint] ;
        }
        return ret;
    }
    //This is not a leaf node
    if (_hasFGTCoefficient) {
        if (_hasFGTCoefficient) {
            double targetError = calculateTargetError(p);
            if (targetError * _error < _rootTree->_epsilon) {
                calculateTargetCoefficients(p);
                ret = 0.;
                for (int iTerm = 0; iTerm < _totalNTerms; iTerm++){
                    ret += _coefficients[iTerm] * targetCoefficients[iTerm];
                }
                ret *= targetError*targetError;
                return ret;
            }
        }
    }
    if (_leftChild != NULL)
        leftValue = _leftChild->fgtEvaluate(p);
    if (_rightChild != NULL)
        rightValue = _rightChild->fgtEvaluate(p);
    ret = leftValue + rightValue;
    return ret;
}

//previously we only calculate the h(y), now we need to calculate h_0--->beta(y) values
void iBoxtree::iBoxtreeNode::fgtEvaluate(double* p,int* beta,int totalNBetas,double* value) {
    //check whether this node can be directly eliminated
    double minimumDistance=0.;
    for(int iDimension=0;iDimension<NDIMENSION;iDimension++){
        double r=0.;
        if(p[iDimension]<_min[iDimension])
            r=_min[iDimension]-p[iDimension];
        else if(p[iDimension]>_max[iDimension])
            r=p[iDimension]-_max[iDimension];
        else
            r=0.;
        minimumDistance+=r*r;
    }
    if(minimumDistance>_rootTree->_cutOffDistance)
        return 0.;
    //Terminal Node
    if (_leftChild == NULL && _rightChild == NULL) {
        if (_hasFGTCoefficient) {
            double targetError = calculateTargetError(p);
            if (targetError * _error < _rootTree->_epsilon) {
                calculateTargetCoefficients(p,beta);
                for (int iTerm = 0; iTerm < totalNBetas; iTerm++) {
                    double ret = 0.;
                    for (int jTerm = 0; jTerm < _totalNTerms; jTerm++) {
                        ret += _coefficients[jTerm] * targetCoefficients[iTerm+jTerm];
                    }
                    ret *= targetError*targetError;
                    value[iTerm]+=ret;
                }
            }
        }
        //we have to directly calculate the summation
        for (int iPoint = _leftPoints; iPoint <= _rightPoints; iPoint++) {
            double temp = 0.;
            for (int iDimension = 0; iDimension < NDIMENSION; iDimension++) {
                double r = (p[iDimension] - _rootTree->_data[iPoint * NDIMENSION + iDimension]) * _rootTree->_invBandwidth;
                temp += r*r;
            }
            temp *= (-0.5);
            temp = exp(temp) * _rootTree->_weight[iPoint];
            int c1 = 1, c2, c3 = 1;
            targetCoefficients[0] = 1.;
            for (int i1 = 0; i1 < NDIMENSION; i1++) {
                double r = (p[i1] - _rootTree->_data[iPoint * NDIMENSION + i1]) * _rootTree->_invBandwidth;
                c2 = c1;
                for (int i3 = c2 - c3; i3 < c2; i3++) {
                    targetCoefficients[c1] = 2. * targetCoefficients[i3] * r;
                    c1++;
                }
                for (int i2 = 2; i2 < beta[i1]; i2++) {
                    c2 = c1;
                    for (int i3 = c2 - c3; i3 < c2; i3++) {
                        targetCoefficients[c1] = 2. * (r * targetCoefficients[i3] - (i2 - 1) * targetCoefficients[i3 - c3]);
                        c1++;
                    }
                }
                c3 *= beta[i1];
            }
            for(int iTerm=0;iTerm<totalNBetas;iTerm++)
                value[iTerm]+=targetCoefficients[iTerm]*temp;
        }
    }
    //This is not a leaf node
    if (_hasFGTCoefficient) {
        if (_hasFGTCoefficient) {
            double targetError = calculateTargetError(p);
            if (targetError * _error < _rootTree->_epsilon) {
                calculateTargetCoefficients(p,beta);
                for (int iTerm = 0; iTerm < totalNBetas; iTerm++) {
                    double ret = 0.;
                    for (int jTerm = 0; jTerm < _totalNTerms; jTerm++) {
                        ret += _coefficients[jTerm] * targetCoefficients[iTerm+jTerm];
                    }
                    ret *= targetError*targetError;
                    value[iTerm]+=ret;
                }
            }
        }
    }
    if (_leftChild != NULL)
        _leftChild->fgtEvaluate(p,beta,value);
    if (_rightChild != NULL)
        _rightChild->fgtEvaluate(p,beta,value);
}

double iBoxtree::iBoxtreeNode::calculateTargetError(double* p) {
    double ret = 1.;
    double r = 0.;
    for (int iDimension = 0; iDimension < NDIMENSION; iDimension++){
        double t=(p[iDimension] - _mean[iDimension])* _rootTree->_invBandwidth;
        r += t*t;
    }
    ret = exp(-0.25 * r);
    return ret;
}

double iBoxtree::hermitePolynomial(double x, int n) {
    double h0 = 1.;
    double ret = 0.;
    if (n == 0)
        return h0;
    double h1 = 2 * x;
    if (n == 1)
        return h1;
    if (n > 1) {
        for (int i1 = 2; i1 <= n; i1++) {
            ret = 2 * (x * h1 - (i1 - 1) * h0);
            h0 = h1;
            h1 = ret;
        }
    }
    return ret;
}

//calculate the hermite coefficients upto _nTerms+nTerms orders
void iBoxtree::iBoxtreeNode::calculateTargetCoefficients(double* p,int* nTerms) {
    int c1 = 1, c2, c3 = 1;
    targetCoefficients[0] = 1.;
    for (int i1 = 0; i1 < NDIMENSION; i1++) {
        double r = (p[i1] - _mean[i1]) * INVERSESQRT2 * _rootTree->_invBandwidth;
        c2 = c1;
        for (int i3 = c2 - c3; i3 < c2; i3++) {
            targetCoefficients[c1]  = 2. * targetCoefficients[i3]*r ;
            c1++;
        }
        for (int i2 = 2; i2 < _nTerms[i1]+nTerms[i1]; i2++) {
            c2 = c1;
            for (int i3 = c2 - c3; i3 < c2; i3++) {
                targetCoefficients[c1] = 2. * (r * targetCoefficients[i3] -  (i2 - 1) * targetCoefficients[i3 - c3]);
                c1++;
            }
        }
        c3 *= (_nTerms[i1]+nTerms[i1]);
    }
//    cout<<p[0]<<", "<<_mean[0]<<endl;
//    for(int iTerm=0;iTerm<_totalNTerms;iTerm++)
//        cout<<_targetCoefficients[iTerm]<<" ";
//    cout<<endl;
//    exit(0);
}

void iBoxtree::iBoxtreeNode::calculateTargetCoefficients(double* p) {
    int c1 = 1, c2, c3 = 1;
    targetCoefficients[0] = 1.;
    for (int i1 = 0; i1 < NDIMENSION; i1++) {
        double r = (p[i1] - _mean[i1]) * INVERSESQRT2 * _rootTree->_invBandwidth;
        c2 = c1;
        for (int i3 = c2 - c3; i3 < c2; i3++) {
            targetCoefficients[c1]  = 2. * targetCoefficients[i3]*r ;
            c1++;
        }
        for (int i2 = 2; i2 < _nTerms[i1]; i2++) {
            c2 = c1;
            for (int i3 = c2 - c3; i3 < c2; i3++) {
                targetCoefficients[c1] = 2. * (r * targetCoefficients[i3] -  (i2 - 1) * targetCoefficients[i3 - c3]);
                c1++;
            }
        }
        c3 *= _nTerms[i1];
    }
//    cout<<p[0]<<", "<<_mean[0]<<endl;
//    for(int iTerm=0;iTerm<_totalNTerms;iTerm++)
//        cout<<_targetCoefficients[iTerm]<<" ";
//    cout<<endl;
//    exit(0);
}
//maximum contribution of this node to the final value
//bandwidth will be chosen to be the maximum
//distance will be chose to be minimum

void iBoxtree::swap(int i, int j) {
    double temp;
    for (int iDimension = 0; iDimension < NDIMENSION; iDimension++) {
        temp = _data[i*NDIMENSION+iDimension];
        _data[i*NDIMENSION+iDimension] = _data[j*NDIMENSION+iDimension];
        _data[j*NDIMENSION+iDimension] = temp;
    }
    temp = _weight[i];
    _weight[i] = _weight[j];
    _weight[j] = temp;
}

void iBoxtree::select(int dimension, int position, int low, int high) {
    int m, r;
    while (low < high) {
        r = (low + high) / 2;
        swap(r, low);
        m = low;
        for (int iPoint = low + 1; iPoint <= high; iPoint++) {
            if (_data[iPoint*NDIMENSION+dimension] < _data[low*NDIMENSION+dimension]) {
                m++;
                swap(m, iPoint);
            }
        }
        swap(low, m);
        if (m <= position) low = m + 1;
        if (m >= position) high = m - 1;
    }
}

void iBoxtree::iBoxtreeNode::buildFGTCoefficient() {
    if (_leftChild)
        _leftChild->buildFGTCoefficient();
    if (_rightChild)
        _rightChild->buildFGTCoefficient();
    if (!_hasFGTCoefficient) return;
    double* tempCoefficients = new double[_totalNTerms];
    if(_coefficients)
        delete[] _coefficients;
    _coefficients            = new double[_totalNTerms];
    memset(_coefficients,0,sizeof(double)*_totalNTerms);
    for (int iPoints = _leftPoints; iPoints <= _rightPoints; iPoints++) {
        int c1 = 1, c2, c3 = 1;
        tempCoefficients[0] = 1.;
        _coefficients[0] += _rootTree->_weight[iPoints];
        for (int i1 = 0; i1 < NDIMENSION; i1++) {
            double r = (_rootTree->_data[iPoints * NDIMENSION + i1] - _mean[i1]) * INVERSESQRT2 * _rootTree-> _invBandwidth;
            for (int i2 = 1; i2 < _nTerms[i1]; i2++) {
                c2 = c1;
                for (int i3 = c2 - c3; i3 < c2; i3++) {
                    tempCoefficients[c1]  =  r * tempCoefficients[i3] / i2;
                    _coefficients[c1]    += (tempCoefficients[c1] * _rootTree->_weight[iPoints]);
                    c1++;
                }
            }
            c3 *= _nTerms[i1];
        }
    }
    delete[] tempCoefficients;
}

void iBoxtree::iBoxtreeNode::buildFGTCoefficient2() {
    if (!_ableBeenApproxmiated) return;
    //calculate the target coefficients
    if(_targetCoefficients)
        delete[] _targetCoefficients;
    _targetCoefficients=new double[_totalNTerms2];
    _rootTree->_root->fgtEvaluate(_mean,_nTerms2,_targetCoefficients);
}
void iBoxtree::buildBox(int low, int high, iBoxtreeNode* root) {
    root->_leftPoints  = low;
    root->_rightPoints = high;
    root->calculateStatistics();
    //check if we need to calculate the fgt coefficients of this node
    double maxDiameter=0;
    int    iMaxDiameter=0;
    for(int iDimension=0;iDimension<NDIMENSION;iDimension++){
        if(root->_diameter[iDimension]>maxDiameter){
            maxDiameter=root->_diameter[iDimension];
            iMaxDiameter=iDimension;
        }
    }
    if(maxDiameter>=_maxX/_invBandwidth){
        root->_error=1.;
        root->_hasFGTCoefficient=false;
    }
    else{
        root->_hasFGTCoefficient=true;
        root->_totalNTerms      =1   ;
        for(int iDimension=0;iDimension<NDIMENSION;iDimension++){
                int ix=int(root->_diameter[iDimension]* _invBandwidth/0.01);
                if(effectiveBins[ix]>_nMaxTerms)
                    root->_nTerms[iDimension]= _nMaxTerms;
                else
                    root->_nTerms[iDimension]= effectiveBins[ix];
                root->_totalNTerms*=root->_nTerms[iDimension]  ;
        }
        double error1 = 1.;
        double error2 = 1.;
        for (int iDimension = 0; iDimension < NDIMENSION; iDimension++) {
            double s     = root->_diameter[iDimension] * _invBandwidth;
            int rlow     = int(ceil(s / _xH));
            int rIndex   = rlow;
            int shift    = (root->_nTerms[iDimension]-_minP)*_maxBeta*_nX;
            error1  *= paramB[rIndex]           ;
            error2  *= paramBP[shift+rIndex]    ;
        }
        root->_error=error1-error2;
    }
    //check if the node needs to be further divided
    //if(root->_error>_epsilon&&root->nPoints()>_nMinimumBoxSize){
    if(root->_error>_epsilon&&root->_totalNTerms<2*root->nPoints()*(NDIMENSION+10)){
        int middle = int(0.5 * (low + high));
        select(iMaxDiameter, middle, low, high);
        root->_cut=0.5*(_data[middle*NDIMENSION+iMaxDiameter]+_data[(middle+1)*NDIMENSION+iMaxDiameter]);
        root->_cutDimension=iMaxDiameter;
        iBoxtreeNode* leftChild = new iBoxtreeNode(this);
        iBoxtreeNode* rightChild = new iBoxtreeNode(this);
        buildBox(low, middle, leftChild);
        leftChild->_parent=root;
        root->_leftChild = leftChild;
        buildBox(middle + 1, high, rightChild);
        rightChild->_parent=root;
        root->_rightChild = rightChild;
    }
    else{
        root->_leftChild=NULL   ;
        root->_rightChild=NULL  ;
    }
    
//    if(!root->_leftChild&&!root->_rightChild){
////        localPCA(root);
//        if(root->_error>1.){
//            cout<<"============================================="<<endl;
//            cout<<"Error "<<root->_error<<" : "<<endl;
//            for(int iDimension=0;iDimension<NDIMENSION;iDimension++){
//                cout<<root->_diameter[iDimension]* _invBandwidth<<" ";
//            }
//            cout<<endl;
//            cout<<"Total NTerms= "<<root->_totalNTerms<<": ";
//            for(int iDimension=0;iDimension<NDIMENSION;iDimension++){
//                cout<<root->_nTerms[iDimension]<<" ";
//            }
//            cout<<endl;
//            for (int iDimension = 0; iDimension < NDIMENSION; iDimension++) {
//                double s = root->_diameter[iDimension] * _invBandwidth;
//                int rlow = int(ceil(s / _xH));
//                int rIndex = _beta[iDimension] * _nX + rlow;
//                int shift = (root->_nTerms[iDimension] - _minP) * _maxBeta*_nX;
//                cout<<"["<<paramB[rIndex]<<", "<<paramBP[shift + rIndex]<<"] ";
//            }
//            cout<<endl;
//        }
//    }
}

double iBoxtree::directEvaluate(double* p){
    double ret=0.       ;
    double Hy =1.       ;
    double h=1.         ;
    for(int iDimension=0;iDimension<NDIMENSION;iDimension++){
        h*=(INVERSESQRT2PI*_invBandwidth);
    }
    for(int iPoint=0;iPoint<_nPoints;iPoint++){
        double r=0.;
        for(int iDimension=0;iDimension<NDIMENSION;iDimension++){
            double s=(p[iDimension]-_data[iPoint*NDIMENSION+iDimension])*_invBandwidth;
            r+=s*s;
        }
        r*=-0.5;
        ret+=exp(r)*Hy*_weight[iPoint];
    }
    return h*ret;
}
double iBoxtree::fgtEvaluate(double* p){
    double ret   =_root->fgtEvaluate(p) ;
    ret         *=_scaleFactor          ;
    return ret;
} 
void iBoxtree::iBoxtreeNode::calculateTargetCoeffcients2(double* p) {
    int c1 = 1, c2, c3 = 1;
    targetCoefficients[0] = 1.;
    for (int i1 = 0; i1 < NDIMENSION; i1++) {
        double r = (p[i1] - _mean[i1]) * INVERSESQRT2 * _rootTree-> _invBandwidth;
        for (int i2 = 1; i2 < _nTerms2[i1]; i2++) {
            c2 = c1;
            for (int i3 = c2 - c3; i3 < c2; i3++) {
                targetCoefficients[c1] += r * targetCoefficients[i3] / i2;
                c1++;
            }
        }
        c3 *= _nTerms2[i1];
    }
}
double iBoxtree::fgtEvaluate2(double* p){
    double ret=0;
    iBoxtreeNode* node=_root;
    while(node){
        if(node->_ableBeenApproxmiated)
            break;
        else if(p[node->_cutDimension]<node->_cut)
            node=node->_leftChild;
        else 
            node=node->_rightChild;
    }
    if(node){
        node->calculateTargetCoeffcients2(p);
        for(int iTerms=0;iTerms<node->_nTerms2;iTerms++){
            ret+=node->_targetCoefficients[iTerms]*targetCoefficients[iTerms];
        }
    }
    else
        ret=fgtEvaluate(p);
    return ret;
}