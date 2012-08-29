//
//  SimpleTest.h
//  LilOpt
//
//  Created by Peter Boyer on 7/29/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef LilOpt_SimpleTest_h
#define LilOpt_SimpleTest_h

#include "Eigen/Dense"
#include "LilOpt"
using namespace LilOpt;
using namespace Eigen;

// Test numeric diff
namespace SimpleTest {
    
    class SimpleFunction : public LilOpt::IErrorFunctionDiff <double, 1, 1, 1> {
        
    public:
        
        virtual bool Evaluate(  const Matrix<double, 1, 1>&    points, 
                                const Matrix<double, 1, 1>&    parms, 
                                Matrix<double, 1, 1>&          residuals,
                                Matrix<double, 1, 1>&          jacobian ) const
        {
            
            residuals = Matrix<double, 1, 1>();
            residuals(0) = 10 - parms(0);
            
            jacobian = Matrix<double, 1, 1>();
            jacobian(0) = -1;
            
            return true;
            
        }
        
        virtual bool Evaluate(  const Matrix<double, 1, 1>&    points, 
                                const Matrix<double, 1, 1>&    parms, 
                                Matrix<double, 1, 1>&          residuals ) const
        {
            
            residuals = Matrix<double, 1, 1>();
            residuals(0) = 10 - parms(0);
            
            return true;
            
        }
        
    };
    
    void test()
    {
        Matrix<double, 1, 1> pts;
            
        IErrorFunctionDiff<double, 1, 1, 1>* errorFunc = new SimpleFunction();
        
        Solver::Options<double> options;
        
        Matrix<double, 1, 1> initParams;
        
        Solver::LevenbergMarquardt<double, 1, 1, 1> lmSolver( options, initParams, errorFunc, pts );
        
        lmSolver.Iterate();
            
    }

}

#endif
