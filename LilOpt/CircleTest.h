//
//  CircleTest.h
//  LilOpt
//
//  Created by Peter Boyer on 7/25/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef LilOpt_CircleTest_h
#define LilOpt_CircleTest_h

#include <time.h>
#include "Eigen/Dense"
#include "IErrorFunction.h"
#include "LevenbergMarquardt.h"
#include "ErrorFunctionNumericDiff.h"
#include "Options.h"

using namespace Eigen;

namespace CircleTest 
{
    
    template<typename _Scalar, unsigned int _NumResiduals>
    Matrix<_Scalar, _NumResiduals, 2> get2DCircleWithNoise( _Scalar noiseMax, 
                                                            _Scalar radius, 
                                                            Matrix<_Scalar,2,1> center)
    
    {
        
        Matrix<_Scalar, _NumResiduals, 2> points;
        
        for (unsigned i = 0; i < _NumResiduals; i++) 
        {
            
            _Scalar f = (_Scalar)rand() / RAND_MAX;
            _Scalar theta = f * 2 * M_PI;
            
            _Scalar f2 = (_Scalar)rand() / RAND_MAX;
            _Scalar noise = f2 * noiseMax - noiseMax / 2;
            
            _Scalar thisRadius = radius + noise;
            
            _Scalar x = center[0] + thisRadius * cos( theta );
            _Scalar y = center[1] + thisRadius * sin( theta );
            
            points.row(i) << x, y;
            
        }
        
        return points;
        
    }
    
    template<typename _Scalar, unsigned int _NumResiduals>
    Matrix<_Scalar, _NumResiduals, 2> get2DRegularCircleWithNoise( _Scalar noiseMax, 
                                                                   _Scalar radius, 
                                                                   Matrix<_Scalar,2,1> center)
    
    {
        
        Matrix<_Scalar, _NumResiduals, 2> points;
        
        for (unsigned i = 0; i < _NumResiduals; i++) 
        {
            
            _Scalar f = (_Scalar)i / (_Scalar)_NumResiduals;
            _Scalar theta = f * 2 * M_PI;
            
            _Scalar f2 = (_Scalar)rand() / RAND_MAX;
            _Scalar noise = f2 * noiseMax - noiseMax / 2;
            
            _Scalar thisRadius = radius + noise;
            
            _Scalar x = center[0] + thisRadius * cos( theta );
            _Scalar y = center[1] + thisRadius * sin( theta );
            
            points.row(i) << x, y;
            
        }
        
        return points;
        
    }
    
    
    
    // _Scalar is templated 
    // _NumParams is 3
    // _NumResiduals is templated
    // _Dimension 2
    
    template<typename _Scalar, unsigned int _NumResiduals >
    class CircleFunction : public LilOpt::IErrorFunction <_Scalar, _NumResiduals, 3, 2> {
        
    public:
        
        virtual bool Evaluate(  const Matrix<_Scalar, _NumResiduals, 2>&    points, 
                                const Matrix<_Scalar, 3, 1>&                parms, 
                                Matrix<_Scalar, _NumResiduals, 1>&          residuals ) const
        {
            
            // go through all the points
            // given the parms, determine the residuals
            // parms are X, Y, R
            
            residuals = Matrix<_Scalar, _NumResiduals, 1>();
            
            Matrix<_Scalar, 1, 2> center = parms.block(0,0,2,1).transpose();
 
            _Scalar radius = parms(2);
            
            for (unsigned i = 0; i < _NumResiduals; i++)
            {
                
                Matrix<_Scalar, 1, 2> point = points.block(i,0,1,2);
                _Scalar dist = (point - center).norm();
                
                if (dist < radius) {
                    residuals[i] = radius - dist;
                } else {
                    residuals[i] = dist - radius;
                }
                
            }
            
            return true;
            
        }
        
    };
    
    void test ()
    {
        // Set up the circle params
        const unsigned int numResids = 20;
        double noiseMax = 0.0;
        double radius = 15.0;
        Matrix<double, 2, 1> center(51,35);
        
        // Generate the circle, given the params
        Matrix<double, numResids, 2> pts = CircleTest::get2DRegularCircleWithNoise<double, numResids>(noiseMax, radius, center);

        // Instantiate the circle function
        LilOpt::IErrorFunction<double, numResids, 3, 2>* circleFunc = new CircleFunction<double, numResids>();
        
        // Create the default options to be given to our solver
        LilOpt::Solver::Options<double> options;
     
        // Create the error function with numeric differentiation for the jacobians
        LilOpt::IErrorFunctionDiff<double, numResids, 3, 2>* funcDiff = new LilOpt::ErrorFunctionNumericDiff<double, numResids, 3, 2>(circleFunc);
        
        // provide the initial params, which we will start from
        Matrix<double, 3, 1> initParams(-20,0,10.0);

        // Create the solver and run it
        LilOpt::Solver::LevenbergMarquardt<double, numResids, 3, 2> lmSolver( options, initParams, funcDiff, pts );
        lmSolver.Minimize();
        
    }
    
    
}


#endif
