LilOpt 0.1.0
======

### Description

LilOpt is a C++, header-based non-linear least squares template library.  Think of it as a light-weight, general multi-dimensional optimization library.  LilOpt provides Levenberg-Marquardt and Gauss-Newton solvers - both industry standard.  


### Current status

LilOpt is undergoing major changes.  Expect the API to change rapidly.


### Future work

Sparse matrices - Once I add support for sparse matrices it will also be much faster for large problems (hundreds or thousands of variables).  Currently, LilOpt uses dense matrices and solving algorithms, which take up prohibitive amounts of memory. 

Multi-threaded solving - LilOpt uses a single threaded Householder QR algorithm for solving the matrices developed on every iteration of the LM and GN algorithm.  I will eventually allow choosing between various algorithms, some of which will provide multi-threading support.  
 
 
### Usage

I've aimed to make it as developer friendly as possible.  It's platform agnostic, although I've included an xcodeproj for anyone who'd like to use it.

LilOpt is the smallest NLLS library I've seen and requires no linking to binaries (as it is header-based).  It's only dependency is the Eigen matrix library, also an open source header template library.  Making it a part of a project is as simple as:

\#include \<LilOpt\>

It is safe to use multiple solvers in separate threads but every solver is in itself single-threaded.


### Examples

There is an example of its usage in CircleTest.h.  In this example, it finds the best fit circle to a set of points (in the least squares sense, of course). 





