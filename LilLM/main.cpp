//
//  main.cpp
//  LilLM
//
//  Created by Peter Boyer on 7/3/12.
//  Copyright (c) 2012 All rights reserved.
//
#include "Eigen/Dense"
#include <iostream>

using namespace Eigen;

// PCA 
// Linear least squares
// Gauss Newton optimization
// Levenberg Marquardt Minimization

namespace LilOpt {

    void Assert(bool predicate, std::string message) {
        
        if (!predicate) {
            
            std::cerr << "LILOPT ERROR" << std::endl;
            std::cerr << message << std::endl;
            exit(1);
            
        }
        
    }
    
    
    // Mostly a warmup.  This class is a subset of the functionality of Eigen.
    
    template<typename _Scalar, unsigned int _Dimension>
    class Linear {
        
    public:
        
        Linear( Matrix<_Scalar, Eigen::Dynamic, _Dimension> dataPoints ):
                _DataPoints(dataPoints)    
        { }
        
        void Regression(    Matrix<_Scalar, Eigen::Dynamic, 1>& b, 
                            Matrix<_Scalar, _Dimension, 1>& solVector ) 
        {

            // Solve
            // A^T * A * x = A^T * b
            Matrix<_Scalar, Eigen::Dynamic, _Dimension>& A = _DataPoints;
            Matrix<_Scalar, _Dimension, Eigen::Dynamic> At = A.transpose();
            
            Matrix<_Scalar, _Dimension, 1> rhs = At*b;
            Matrix<_Scalar, _Dimension, _Dimension> lhs = At*A;
    
            solVector = lhs.fullPivHouseholderQr().solve( rhs );
            
        }
        
    private:
        
        Eigen::Matrix<_Scalar, Eigen::Dynamic, _Dimension> _DataPoints;
        
    };
    

}    
    






class GaussNewton {
    
    GaussNewton( ) {
        
        
        
    }
    
};

int DoItBaby(int n) {
    int w = n*n;
    return w;
}

int main(int argc, const char * argv[])
{

    Matrix<float, 3, 2> vals;
    vals << 1, -2, 
            1, 0, 
            1, 2;
    
    LilOpt::Linear<float, 2> solver( vals );
    
    Matrix<float, Eigen::Dynamic, 1> b;
    
    b << 1, 
    2, 4;
    
    Matrix<float, 2, 1> xHat;
    solver.Regression(b, xHat);
    
    std::cout << xHat << std::endl;
    
    return 0;
}


// these are nice but should come later
template<typename _Scalar, int _NumParams>
Matrix<_Scalar, 1, _NumParams> Jacobian1d( Matrix<_Scalar, _NumParams, 1> parms) {
    
    
    
    return Matrix<_Scalar, 1, _NumParams>();
    
}

template<typename _Scalar, int _NumParams>
_Scalar ErrorFunc1d( Matrix<_Scalar, _NumParams, 1> parms ) {
    
    _Scalar f = 0;
    
    return f;
    
}