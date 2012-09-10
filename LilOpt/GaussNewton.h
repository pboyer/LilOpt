//
//  GaussNewton.h
//  LilOpt
//
//  Created by Peter Boyer on 7/9/12.
//  Copyright (c) 2012 All rights reserved.
//

#ifndef LilOpt_GaussNewton_h
#define LilOpt_GaussNewton_h

#include "Options.h"
#include "Eigen/Dense"

using namespace Eigen;

namespace LilOpt {
    namespace Solver {
        
        template< typename _Scalar, unsigned int _NumResiduals, unsigned int _NumParams, unsigned int _Dimension >
        class GaussNewton {
            
        public:
            
        // CONSTRUCTOR /////////////////////////
            
            GaussNewton(    Options<_Scalar> options,
                            Matrix<_Scalar, _NumParams, 1>& initialParams,
                            IErrorFunctionDiff<_Scalar, _NumResiduals, _NumParams, _Dimension>* function, 
                            Matrix<_Scalar, Eigen::Dynamic, _Dimension>& dataPoints)
                            :
                            _SumResidual(-1),
                            _LastSumResidual(-1),
                            _Options(options), 
                            _DataPoints(dataPoints), 
                            _Function(function), 
                            _CurrentParams(initialParams)  
            
            {}
            
        // MEMBER METHODS //////////////////////
            
            bool Iterate() {
                
                // Evaluate residuals and jacobian at _Current Params
                Matrix<_Scalar, _NumResiduals, _NumParams> J;
                Matrix<_Scalar, _NumResiduals, 1>          Ep0;
                _Function->Evaluate(  _DataPoints, _CurrentParams, Ep0, J );
                UpdateSumResiduals( Ep0.sum() );
                
                // Solve.  
                // J^T * J ( P1 - P0 ) = J^T * Ep0
                Matrix<_Scalar, _NumParams, _NumResiduals> Jt = J.transpose();
                Matrix<_Scalar, _NumParams, 1> P1mP0 = ( Jt * J ).fullPivHouseholderQr().solve( - Jt * Ep0 );  
                
                _CurrentParams = P1mP0 + _CurrentParams;
                
                return true;
            }
            
            bool Minimize() 
            {
                
                // Iterate upto MaxIterations, terminate if difference between 
                // last and new sumResidual is less than value
                unsigned i = 0;
                bool iterationSuccess = true;
                
                do {
                    
                    iterationSuccess = Iterate();
                    
                } while (i < _Options.MaxIterations 
                         && iterationSuccess 
                         && abs( _LastSumResidual - _SumResidual) > _Options.Tolerance );
                
                return true;
                
            }
            
        // PRIVATE METHODS ///////////////////
            
        private:
            
            void UpdateSumResiduals( _Scalar newSumResidual )
            {
                _LastSumResidual = _SumResidual;
                _SumResidual = newSumResidual;
            }
            
        // PUBLIC FIELDS //////////////////////
            
            Options<_Scalar> _Options;
            
        // PRIVATE FIELDS /////////////////////
            
        private:
            
            _Scalar         _SumResidual;
            _Scalar         _LastSumResidual;
            IErrorFunctionDiff<_Scalar, _NumResiduals, _NumParams, _Dimension>* _Function;
            Eigen::Matrix<_Scalar, _NumResiduals, _Dimension> _DataPoints;
            Eigen::Matrix<_Scalar, _NumParams, 1> _CurrentParams;
            
        };
    }
}

#endif
