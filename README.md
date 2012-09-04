LilOpt
======

LilOpt is a C++, header-based non-linear least squares template library.  Think of it as a light-weight, general multi-dimensional optimization library.  

LilOpt provides Levenberg-Marquardt and Gauss-Newton solvers - both industry standard.  Once I add support for sparse matrices it will also be much faster for large problems (hundreds or thousands of variables).  Currently, it uses dense matrices and solving algorithms, which take up prohibitive amounts of memory.  

There is an example of its usage in CircleTest.h.  In this example, it finds the best fit circle to a set of points (in the least squares sense, of course).  

I've aimed to make it as developer friendly as possible.  It's platform agnostic, although I've used xcode here (should migrate to cmake).  It's the smallest NLLS library I've seen and requires no linking to binaries (as it is header-based).  It's only dependency is the Eigen matrix library, also an open source header template library.  Making it a part of a project is as simple as:

"#include <LilOpt>"

