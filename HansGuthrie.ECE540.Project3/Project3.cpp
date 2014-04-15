#include "matrix.hpp"
#include <stdio.h>

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#ifndef MatrixRead_h
#define MatrixRead_h 0


matrix ReadBinaryMatrix( char *Name )
{
	int Rows, Columns;
	FILE *fin;
	matrix output;
	fin = fopen( Name, "rb" );

	//Gets the rows/cols from the beginning of the file
	fread( &Rows, sizeof( int ), 1, fin );
	fread( &Columns, sizeof( int ), 1, fin );

	output = matrix( Rows, Columns ); //creates  the output matrix

	if ( fin ) //if successfully opened file for reading
	{
		fread( output.AsPointer( ), sizeof( double ), Rows*Columns, fin ); //Read through the file for elements the size of a double
		fclose( fin );
	}
	return output;
} //end ReadBinaryMatrix

#endif

//LinearRegressionAndCC.pdf
// Regression analysis using two-pass Approach, return vector of parameters p.
// y(i) = p(0) * x(i) + p(1) .
matrix FaultTolerantRegression( matrix x, matrix y )
{
	int k;
	double x_bar = 0, y_bar = 0,
		Numerator = 0, Denominator = 0;
	// Test that the two vectors are the same length.
	if ( x.high( ) != y.high( ) )
	{
		matrix z;
		return z;
	} // End of invalid input test.
	matrix z( 2 ); // Create output matrix.
	// Loop through to compute mean of x and y.
	for ( k = 0; k < x.high( ); k++ )
	{
		x_bar += x( k );
		y_bar += y( k );
	} // End of loop to compute means.
	// Scale sums to be means.
	x_bar /= (double)x.high( );
	y_bar /= (double)x.high( );
	// Loop computing next step in computing a
	for ( k = 0; k < x.high( ); k++ )
	{
		Numerator += ( x( k ) - x_bar )*( y( k ) - y_bar );
		Denominator += ( x( k ) - x_bar )*( x( k ) - x_bar );
	} // End of loop to compute a
	// Compute offset,
	z( 0 ) = Numerator / Denominator;
	z( 1 ) = y_bar - z( 0 ) * x_bar;
	return z;
} // End of FaultTolerantRegression


// support for **CorrelationCoefficient**
double SumOfArraysMinusMeansMultiplied( double *x, double *y, int length,
	double xMean, double yMean )
{
	double sum = ( *x - xMean ) * ( *y - yMean );
	while ( --length )
	{
		x++; y++;
		sum += ( *x - xMean ) * ( *y - yMean );
	} // end of while loop.
	return sum;
}// end of SumOfArraysMinusMeansMultiplied

// Code for computing the Correlation Coefficient
double CorrelationCoefficient( matrix x, matrix y ) // input vectors.
{
	double cc, xMean, yMean;
	double xSigma, ySigma;
	if ( x.high( ) != y.high( ) )
	{
		return NAN;
	}
	xMean = ComputeMean( x.AsPointer( ), x.high( ) );
	yMean = ComputeMean( y.AsPointer( ), y.high( ) );
	xSigma = ComputeStdev( x.AsPointer( ), x.high( ), xMean );
	ySigma = ComputeStdev( y.AsPointer( ), y.high( ), yMean );
	cc = SumOfArraysMinusMeansMultiplied( x.AsPointer( ), y.AsPointer( ), y.high( ),
		xMean, yMean ) / (double)( x.high( ) - 1 );
	return cc / ( xSigma * ySigma );
} // end of CorrelationCoefficient

void main( )
{
	matrix A;
	A = ReadBinaryMatrix( "StatsData.mtx" );

}