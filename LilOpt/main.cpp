#include "Eigen/Dense"
#include "LilOpt"
#include <iostream>
#include <utility>

using namespace Eigen;

//
// Lilopt:  A C++ header library for non-linear least squares optimization
//

void LinearRegressionTestRowsNotMatched() {
    
    Matrix<float, 4, 2> vals;
    
    vals << 1, -2, 
    1,  0, 
    1,  2,
    4,  5;
    
    LilOpt::Solver::Linear<float, 2> solver( vals );
    
    Matrix<float, Eigen::Dynamic, 1> b(3);
    
    b <<    1, 
    2, 
    4;
    
    Matrix<float, 2, 1> xHat;
    solver.Regress( b, xHat);
    
}

void LinearRegressionTestRowsMatched() {
    
    Matrix<float, 4, 2> vals;
    
    vals << 1, -2, 
    1,  0, 
    1,  2;
    
    LilOpt::Solver::Linear<float, 2> solver( vals );
    
    Matrix<float, Eigen::Dynamic, 1> b(3);
    
    b <<    1, 
    2, 
    4;
    
    Matrix<float, 2, 1> xHat;
    solver.Regress( b, xHat);
    
    std::cout << xHat << std::endl;
    
}





int main(int argc, const char * argv[])
{
    
    
    
    
    return 0;
}