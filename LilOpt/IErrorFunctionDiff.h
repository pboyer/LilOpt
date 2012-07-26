//
//  ErrorFunctionDiff.h
//  LilOpt
//
//  Created by Peter Boyer on 7/9/12.
//  Copyright (c) 2012 All rights reserved.
//

#ifndef LilLM_IErrorFunctionDiff_h
#define LilLM_IErrorFunctionDiff_h

#include "Eigen/Dense"

using namespace Eigen;

namespace LilOpt {
    
    template<typename _Scalar, unsigned int _NumResiduals, unsigned int _NumParams, unsigned int _Dimension>
    class IErrorFunctionDiff {
        
    public:
        
        virtual bool Evaluate(  const Matrix<_Scalar, _NumResiduals, _Dimension>&   points, 
                                const Matrix<_Scalar, _NumParams, 1>&               parms, 
                                Matrix<_Scalar, _NumResiduals, 1>&                  residuals, 
                                Matrix<_Scalar, _NumResiduals, _NumParams>&         jacobian ) const = 0;
        
        virtual bool Evaluate(  const Matrix<_Scalar, _NumResiduals, _Dimension>&   points, 
                                const Matrix<_Scalar, _NumParams, 1>&                parms, 
                                Matrix<_Scalar, _NumResiduals, 1>&                  residuals) const = 0;
        
    };
    
}

#endif
