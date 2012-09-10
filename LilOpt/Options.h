//
//  Options.h
//  LilOpt
//
//  Created by Peter Boyer on 7/9/12.
//  Copyright (c) 2012 All rights reserved.
//

#ifndef LilOpt_Options_h
#define LilOpt_Options_h

namespace LilOpt {
    namespace Solver {

        enum SolverType { DENSE_QR };
        
        template<typename _Scalar>
        struct Options {
        
        public:
        
            _Scalar         Tolerance;  
            unsigned int    MaxIterations;
            SolverType      SolverType;
            bool            WriteProgressToStdout;
            _Scalar         LevenbergMarquardtLambda;
            _Scalar         LevenbergMarquardtV;
            unsigned int    MaxSubIterations;
        
            Options():  MaxIterations(100), 
                        SolverType( DENSE_QR ), 
                        WriteProgressToStdout(true),
                        LevenbergMarquardtLambda(100),
                        LevenbergMarquardtV(10),
                        Tolerance(1e-2),
                        MaxSubIterations(10)
            
            {}
            
            
            
        };
    }
}

#endif
