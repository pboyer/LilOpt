//
//  ErrorFunction.h
//  LilOpt
//
//  Created by Peter Boyer on 7/9/12.
//  Copyright (c) 2012 All rights reserved.
//

#ifndef LilLM_IErrorFunction_h
#define LilLM_IErrorFunction_h

#include "Eigen/Dense"

using namespace Eigen;

namespace LilOpt {

    template<typename _Scalar, unsigned int _NumResiduals, unsigned int _NumParams, unsigned int _Dimension>
    class IErrorFunction {
    
    public:
    
        virtual bool Evaluate(  const Matrix<_Scalar, _NumResiduals, _Dimension>&    points, 
                                const Matrix<_Scalar, _NumParams, 1>&                parms, 
                                Matrix<_Scalar, _NumResiduals, 1>&                   residuals ) const = 0;
    
    };
    
}

#endif
