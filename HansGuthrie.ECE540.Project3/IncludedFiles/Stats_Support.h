#ifndef Stats_Support_h
#define Stats_Support_h

#include "Matrix.hpp"

#ifndef NAN
#define NAN_HUGE_ENUF  1e+300	/* _HUGE_ENUF*_HUGE_ENUF must overflow */
#define NAN_INFINITY   ((float)(NAN_HUGE_ENUF * NAN_HUGE_ENUF))
#define NAN  ((float)(NAN_INFINITY * 0.0F))
#endif
// Regression analysis using two-pass Approach, return coefficient in vector p.
//   y(i) = p(0) * x(i) + p(1) ...
matrix FaultTolerantRegression(matrix x, matrix y);

// This function will compute the value of the linear equation defined
// by parms for each value in the array x, placing the result in the return vector.   
matrix ComputeLinear(matrix x, matrix p);
// Code for computing the Correlation Coefficient 
double CorrelationCoefficient(matrix x, matrix y);  // input vectors.

// Function to computer coefficient of determination
double CoefficientOfDetermination(matrix y, matrix z);

#endif
