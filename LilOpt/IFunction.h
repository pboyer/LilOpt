//
//  Function.h
//  LilOpt
//
//  Created by Peter Boyer on 7/9/12.
//  Copyright (c) 2012 All rights reserved.
//

#ifndef LilLM_IFunction_h
#define LilLM_IFunction_h

#include "Eigen/Dense"

using namespace Eigen;

namespace LilOpt {
    
    template<typename _Scalar, unsigned int _NumParams, unsigned int _Dimension>
    class IFunction {
        
    public:
    
        virtual bool Evaluate(  const Matrix<_Scalar, _NumParams, 1>& parms, 
                                Matrix<_Scalar, _Dimension, 1> & value ) const = 0;
    
    };
    
}

#endif
