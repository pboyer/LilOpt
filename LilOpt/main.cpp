#include "Eigen/Dense"
#include "LilOpt"
#include <iostream>
#include <utility>
#include <vector>

#include "CircleTest.h"
#include "SimpleTest.h"

using namespace Eigen;
using namespace std;

//
// Lilopt:  A C++ header library for non-linear least squares optimization
//


void LinearRegressionTestRowsMatched() {
    
    Matrix<float, 3, 2> vals;
    
    vals << 1, -2, 
    1,  0, 
    1,  2;
    
    LilOpt::Solver::Linear<float, 3, 2> solver( vals );
    
    Matrix<float, 3, 1> b;
    
    b <<    1, 
    2, 
    4;
    
    Matrix<float, 2, 1> xHat;
    solver.Regress( b, xHat);
    
    std::cout << xHat << std::endl;
    
}

int main(int argc, const char * argv[])
{
    
    //SimpleTest::test();
    CircleTest::test();
    
    
    return 0;
}




