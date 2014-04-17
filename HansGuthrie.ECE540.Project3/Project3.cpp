#include "matrix.hpp"
#include <stdio.h>
#include "MatrixRead.h"
#include "Stats_Support.h"
#include "MatrixOutputs.hpp"
#include <vector>
#include "Project3.h"

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif


///Computes the Multivariable regression of matrix a, matrix b, and the dependent matrix. 
///Requires input of the size of the matricies as an int
void MultiVariableRegression( matrix a, matrix b, matrix dependent, int size )
{
	matrix In = matrix( size, 3 );

	matrix PIn;
	matrix p, x;
	double CoD, sqrootofCoD, CC, temp;

	//create matrix In from Row 3, Row 4, and a lot of 1's.
	for ( int k = 0; k < size; k++ )
	{
		In( k, 0 ) = a( k ); // row 3
		In( k, 1 ) = b( k ); //row 4
		In( k, 2 ) = 1.0; //set all of last column to 1's
	}

	PIn = MatrixPseudoInverse( In );
	if ( PIn.isValid( ) )
	{
		// compute the pseudo inverse and check for valid operation.
		// Multiply pseudo inverse and dependent (row5), and check for valid operation.
		p = PIn * dependent;
		if ( p.isValid( ) )
		{
			//prints p
			PrintMatrix( p );

			x = In * p; //x matrix is the product of In and p
			if ( x.isValid( ) )
			{
				CoD = CoefficientOfDetermination( x, dependent );
				sqrootofCoD = sqrt( CoD );
				CC = CorrelationCoefficient( x, dependent );
				printf( "Correlation of Determination (CoD) = %lg\n", CoD );
				printf( "Sqrt of CoD = %lg\n", sqrootofCoD );
				printf( "Correlation Coefficient (CC) = %lg\n", CC );

			}// end of valid x = MtrxVector check
		} // end of valid p = MtrxVector check.
	}
}

void ConfidenceInterval( matrix m, int samplesize )
{
	//mean, variance
	double Mean, Stdev, MeanLow, MeanHigh;
	int Length = m.high( );

	Mean = ComputeMean( m.AsPointer( ), Length );
	Stdev = ComputeStdev( m.AsPointer( ), Length, Mean );
	MeanLow = Mean - 1.64*Stdev / sqrt( (double)samplesize );
	MeanHigh = Mean + 1.64*Stdev / sqrt( (double)samplesize );

	printf( "Mean: %lf\nStdev: %lf\n", Mean, Stdev );
	printf( "90%% Confidence of being between %lf and %lf\n\n", MeanLow, MeanHigh );
}


int main( )
{
	matrix InputMatrix; //the main matrix that we're reading in
	matrix row1, row2, row3, row4, row5; //separate matricies for each row
	matrix FaultTol1, FaultTol2, FaultTol3, FaultTol4;
	double CorrelCoeff1, CorrelCoeff2, CorrelCoeff3, CorrelCoeff4;

	//1a)

	InputMatrix = ReadBinaryMatrix( "StatsData.mtx" ); //read in the matrix from file
	//calculate the boundaries for the matricies so we don't have to calculate it multiple times
	int numrows = InputMatrix.high( );
	int numcols = InputMatrix.wide( );

	row1 = matrix( numcols, 1 );
	row2 = matrix( numcols, 1 );
	row3 = matrix( numcols, 1 );
	row4 = matrix( numcols, 1 );
	row5 = matrix( numcols, 1 );

	//matrix is created in this format (rows, cols)
	//(row1 = InputMatrix row 0)
	for ( int i = 0; i < numcols; i++ )
	{
		row1( i ) = InputMatrix( 0, i );
	}
	for ( int i = 0; i < numcols; i++ )
	{
		row2( i ) = InputMatrix( 1, i );
	}
	for ( int i = 0; i < numcols; i++ )
	{
		row3( i ) = InputMatrix( 2, i );
	}
	for ( int i = 0; i < numcols; i++ )
	{
		row4( i ) = InputMatrix( 3, i );
	}
	for ( int i = 0; i < numcols; i++ )
	{
		row5( i ) = InputMatrix( 4, i );
	}

	//Compute Linear Regression parameters (FaultTol)
	FaultTol1 = FaultTolerantRegression( row1, row5 );
	FaultTol2 = FaultTolerantRegression( row2, row5 );
	FaultTol3 = FaultTolerantRegression( row3, row5 );
	FaultTol4 = FaultTolerantRegression( row4, row5 );

	//Compute Correlation Coeffcient 
	CorrelCoeff1 = CorrelationCoefficient( row1, row5 );
	CorrelCoeff2 = CorrelationCoefficient( row2, row5 );
	CorrelCoeff3 = CorrelationCoefficient( row3, row5 );
	CorrelCoeff4 = CorrelationCoefficient( row4, row5 );

	//Print what we have computed

	printf( "Linear Regression parameters\n" );
	PrintMatrix( FaultTol1 );
	printf( "\n" );
	PrintMatrix( FaultTol2 );
	printf( "\n" );
	PrintMatrix( FaultTol3 );
	printf( "\n" );
	PrintMatrix( FaultTol4 );
	printf( "\n" );

	printf( "Correlation Coefficients\n" );
	printf( "CC1: %lg \n", CorrelCoeff1 );
	printf( "CC2: %lg \n", CorrelCoeff2 );
	printf( "CC3: %lg \n", CorrelCoeff3 );
	printf( "CC4: %lg \n", CorrelCoeff4 );

	//1b)

	//Row 3 and Row 4 have the highest CC, so perform a multivariable regression on them
	MultiVariableRegression( row3, row4, row5, numcols );


	//2a)
	printf( "\n" );
	int samplesize = 81;
	printf( "Row 1 stats:\n" );
	ConfidenceInterval( row1, samplesize );
	printf( "Row 2 stats:\n" );
	ConfidenceInterval( row2, samplesize );
	printf( "Row 3 stats:\n" );
	ConfidenceInterval( row3, samplesize );
	printf( "Row 4 stats:\n" );
	ConfidenceInterval( row4, samplesize );

	

	getchar( );

	return 0;
}
