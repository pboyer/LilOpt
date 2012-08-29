//
//  LevenbergMarquardt.h
//  LilOpt
//
//  Created by Peter Boyer on 7/9/12.
//  Copyright (c) 2012 All rights reserved.
//

#ifndef LilOpt_LevenbergMarquardt_h
#define LilOpt_LevenbergMarquardt_h

#include "Options.h"
#include "Eigen/Dense"
#include "IErrorFunctionDiff.h"
#include <iostream>

using namespace Eigen;

namespace LilOpt {
    namespace Solver {
        
        template< typename _Scalar, unsigned int _NumResiduals, unsigned int _NumParams, unsigned int _Dimension >
        class LevenbergMarquardt {
            
        public:
            
        // CONSTRUCTOR /////////////////////////
            
            LevenbergMarquardt( Options<_Scalar> options,
                                Matrix<_Scalar, _NumParams, 1>& initialParams,
                                IErrorFunctionDiff<_Scalar, _NumResiduals, _NumParams, _Dimension>* function, 
                                Matrix<_Scalar, _NumResiduals, _Dimension>& dataPoints):
            
                                _Options(options), 
                                _DataPoints(dataPoints), 
                                _Function(function), 
                                _CurrentParams(initialParams),
                                _SumResidual(-1) // We haven't run an iteration.  TODO: this is a hack...
            {}
            
        // MEMBER METHODS //////////////////////
            
            bool Iterate() 
            {
                
                // Evaluate residuals and jacobian at _CurrentParams
                Matrix<_Scalar, _NumResiduals, _NumParams> J;
                Matrix<_Scalar, _NumResiduals, 1>          Ep0;
                _Function->Evaluate( _DataPoints, _CurrentParams, Ep0, J );
                
                std::cout << "J" << std::endl;
                std::cout << J << std::endl << std::endl;
                std::cout << "Ep0" << std::endl;
                std::cout << Ep0 << std::endl << std::endl;
                
                // Solve for P1
                // (J^T * J + a * diag(J^T * J) )( P1 - P0 ) = J^T * Ep0
                
                Matrix<_Scalar, _NumParams, _NumResiduals> Jt = J.transpose();
                Matrix<_Scalar, _NumParams, _NumParams> JtJ = Jt * J;
                Matrix<_Scalar, _NumParams, _NumParams> JtJdiag = JtJ.diagonal().asDiagonal();
                Matrix<_Scalar, _NumParams, 1> JtEp0 = -Jt * Ep0;
                Matrix<_Scalar, _NumParams, 1> P1mP0 = ( JtJ  + _Options.LevenbergMarquardtLambda * JtJdiag ).fullPivHouseholderQr().solve( JtEp0 );
                Matrix<_Scalar, _NumParams, 1> P1 = P1mP0 + _CurrentParams;
                
                // evaluate the residual at the new parameter position
                
                Matrix<_Scalar, _NumResiduals, 1>          Ep1;
                _Function->Evaluate( _DataPoints, P1, Ep1 );
                
                // Compute residual sums at start and end parameters
                
                _Scalar Ep0Sum = Ep0.sum();
                _Scalar Ep1Sum = Ep1.sum();
                
             //  DEBUG
//                std::cout << "JtJDiag: "<< std::endl << JtJdiag << std::endl;
                std::cout << "Ep0Sum: " << Ep0Sum << std::endl;
                std::cout << "Ep1Sum: " << Ep1Sum << std::endl;
                std::cout << "P1: " << std::endl << P1 << std::endl << std::endl;
                
                _Scalar v = _Options.LevenbergMarquardtV;
                
                unsigned its = 0;
                
                while (Ep1Sum > Ep0Sum && its < _Options.MaxSubIterations ) 
                {
                    
                    its++;
                    
                    // if Ep1Sum > Ep0Sum, we increase a and reevaluate
                    _Options.LevenbergMarquardtLambda *= v;
                    
                    // Resolve for P1 with new lambda
                    P1mP0 = ( JtJ  + _Options.LevenbergMarquardtLambda * JtJdiag ).fullPivHouseholderQr().solve( JtEp0 );
                    P1 = P1mP0 + _CurrentParams;
                    
                    _Function->Evaluate( _DataPoints, P1, Ep1 );
                    Ep1Sum = Ep1.sum();
                    
                //  DEBUG
                    std::cout << "Sub-Iteration: " << its << std::endl;
                    std::cout << "Ep1Sum: " << Ep1Sum << std::endl;
                    std::cout << "P1: " << P1 << std::endl<< std::endl;
                    
                }
                
                _CurrentParams = P1;
                _Options.LevenbergMarquardtLambda /= v;  
                
                return true;
            }
            
            bool Minimize() 
            {
                
                // Iterate upto MaxIterations, terminate if difference between 
                // last and new sumResidual is less than value
                for (unsigned int i = 0; i < _Options.MaxIterations; i++) {
                    bool result = Iterate();
                    // TODO: check with tolerance
                    if (!result)
                        return result;
                }
                
                return true;
            }
            
        // PUBLIC FIELDS //////////////////////
            
            Options<_Scalar> _Options;
            
        // PRIVATE FIELDS /////////////////////
            
        private:
            
            _Scalar _SumResidual;
            IErrorFunctionDiff<_Scalar, _NumResiduals, _NumParams, _Dimension>* _Function;
            Matrix<_Scalar, _NumResiduals, _Dimension> _DataPoints;
            Matrix<_Scalar, _NumParams, 1> _CurrentParams;
            
        };
        
    } // namespace Solver
}

#endif
