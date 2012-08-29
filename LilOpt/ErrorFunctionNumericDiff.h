//
//  ErrorFunctionNumericDiff.h
//  LilOpt
//
//  Created by Peter Boyer on 7/25/12.
//  Copyright (c) 2012 All rights reserved.
//

#ifndef LilLM_ErrorFunctionNumericDiff_h
#define LilLM_ErrorFunctionNumericDiff_h

namespace LilOpt {

    // TODO: Implement forward, backward
    enum NumericDiffType { CENTERED };

    // Implements a function to perform centered difference evaluation of the jacobian
    template<   typename _Scalar, 
                unsigned int _NumResiduals, 
                unsigned int _NumParams, 
                unsigned int _Dimension >
    class ErrorFunctionNumericDiff : public IErrorFunctionDiff< _Scalar, _NumResiduals, _NumParams, _Dimension > {
    
    public: 
        
    // CONSTRUCTOR /////////////////////////
        
        ErrorFunctionNumericDiff( IErrorFunction<_Scalar, _NumResiduals, _NumParams, _Dimension >* errorFunction ) {
            
            _ErrorFunction = errorFunction;
            
        }
        
    // MEMBER METHODS //////////////////////
        
        virtual bool Evaluate(  const Matrix<_Scalar, _NumResiduals, _Dimension>&       points, 
                                const Matrix<_Scalar, _NumParams, 1>&                   parms, 
                                Matrix<_Scalar, _NumResiduals, 1>&                      residuals, 
                                Matrix<_Scalar, _NumResiduals, _NumParams>&             jacobian ) const
        {
            
            const _Scalar stepSize = 1e-10;
            const _Scalar stepSizeX2 = 2 * stepSize;
            
            _ErrorFunction->Evaluate(points, parms, residuals);
            
            jacobian = Matrix<_Scalar, _NumResiduals, _NumParams>();
            
            for (unsigned i = 0; i < _NumParams; i++) {
                
                Matrix<_Scalar, _NumParams, 1> parmsForw = parms;
                Matrix<_Scalar, _NumResiduals, 1> residualForw;
                parmsForw[i] = parmsForw[i] + stepSize;
                
                Matrix<_Scalar, _NumParams, 1> parmsBackw = parms;
                Matrix<_Scalar, _NumResiduals, 1> residualBackw;
                parmsBackw[i] = parmsBackw[i] - stepSize;
                
                // Perform forward evaluation of residuals
                _ErrorFunction->Evaluate(points, parmsForw, residualForw);
                
                // Perform backward evaluation of residuals
                _ErrorFunction->Evaluate(points, parmsBackw, residualBackw);
                
                // Determine centered finite difference and store in gradient
                jacobian.col(i) = (residualForw - residualBackw) / stepSizeX2;
                
            }  
            
            return true;
            
        }
        
        // Compute the residuals at the given parameters, don't compute the jacobian
        virtual bool Evaluate(  const Matrix<_Scalar, _NumResiduals, _Dimension>&       points, 
                                const Matrix<_Scalar, _NumParams, 1>&                   parms, 
                                Matrix<_Scalar, _NumResiduals, 1>&                      residuals ) const
        {
            return _ErrorFunction->Evaluate( points, parms, residuals);
        }
        
        
    // PRIVATE FIELDS /////////////////////
        
    private:
    
        IErrorFunction<_Scalar, _NumResiduals, _NumParams, _Dimension >* _ErrorFunction;
        
    };

}

#endif
