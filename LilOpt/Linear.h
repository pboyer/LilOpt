//
//
//  Linear.h
//  LilOpt
//
//  Created by Peter Boyer on 7/9/12.
//  Copyright (c) 2012 All rights reserved.
//

#ifndef LilOpt_Linear_h
#define LilOpt_Linear_h

#include "Eigen/Dense"
#include "Assert.h"

using namespace Eigen;

namespace LilOpt {
    
    namespace Solver {

        template<typename _Scalar, unsigned int _NumResiduals, unsigned int _Dimension>
        class Linear {
            
        public:
    
            Linear( Matrix<_Scalar, _NumResiduals, _Dimension> dataPoints ):
            _DataPoints(dataPoints)    
            { }
    
            bool Regress(   Matrix<_Scalar, _NumResiduals, 1>& b, 
                            Matrix<_Scalar, _Dimension, 1>& solVector ) 
            {
        
                // Solve
                // A^T * A * x = A^T * b
                Matrix<_Scalar, _NumResiduals, _Dimension>& A = _DataPoints;
                Matrix<_Scalar, _Dimension, _NumResiduals> At = A.transpose();
                Matrix<_Scalar, _Dimension, 1> rhs = At*b;
                Matrix<_Scalar, _Dimension, _Dimension> lhs = At*A;
        
                solVector = lhs.fullPivHouseholderQr().solve( rhs );
        
                return true;
            }
    
        private:
    
            Eigen::Matrix<_Scalar, _NumResiduals, _Dimension> _DataPoints;
    
        };
    
    }
    
}

#endif
