#include "matrix.hpp"
#include <stdio.h>
#include "MatrixRead.h"
#include "Stats_Support.h"
#include "MatrixOutputs.hpp"
#include <vector>
#include <math.h>
#include "Project3.h"
#include "RandomNumbers.h"

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

///<Summary>
///Computes the Multivariable regression of matrix a, matrix b, and the dependent matrix.
///Requires input of the size of the matricies as an int
///</Summary>
void MultiVariableRegression( matrix a, matrix b, matrix dependent, int size )
{
	matrix In = matrix( size, 3 );

	matrix PIn;
	matrix p, x;
	double CoD, sqrootofCoD, CC;

	//create matrix In from Row 3, Row 4, and a lot of 1's.
	for ( int k = 0; k < size; k++ )
	{
		In( k, 0 ) = a( k ); // row 3
		In( k, 1 ) = b( k ); //row 4
		In( k, 2 ) = 1.0; //set all of last column to 1's
	}

	PIn = MatrixPseudoInverse( In ); //calculate the Pseudoinverse
	if ( PIn.isValid( ) )
	{
		// compute the pseudo inverse and check for valid operation.
		// Multiply pseudo inverse and dependent (row5), and check for valid operation.
		p = PIn * dependent;
		if ( p.isValid( ) )
		{
			//prints p
			printf( "Parameters" );
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
} //end MultiVariableRegression

///<Summary>
///Calculates the confidence internal for the m matrix.
///Requires an int representing the number of items in the sample size (for us 81)
///</Summary>
///<param name="m">The matrix</param>
///<param name="samplesize">Num items in bin</param>
void ConfidenceInterval( matrix m, int samplesize )
{
	double Mean, Stdev, MeanLow, MeanHigh;
	int Length = m.high( );

	Mean = ComputeMean( m.AsPointer( ), Length );
	Stdev = ComputeStdev( m.AsPointer( ), Length, Mean );
	MeanLow = Mean - 1.64*Stdev / sqrt( (double)samplesize );
	MeanHigh = Mean + 1.64*Stdev / sqrt( (double)samplesize );

	printf( "Mean: %lf\nStdev: %lf\n", Mean, Stdev );
	printf( "90%% Confidence of being between %lf and %lf\n\n", MeanLow, MeanHigh );
} //end ConfidenceInterval

///<summary>
///Computes the mean and the standard deviation of a pointer to a row array.
///Prints the Mean, Standard Deviation, Number of times SIMs inside CI:
///</summary>
///<param name = "*row">Pointer to the row's data (double array)</param>
void IntervalMeanAndDev( double *row )
{
	const int binsize = 2000; //the number of subsections in the matrix
	double mean[ binsize ], stdev[ binsize ]; //arrays to hold the mean and stdev for what we will compute
	double ComputedMean, ComputedStdev, CIHigh, CILow;
	int NumTimesInsideCI = 0;

	//loop through each bin and compute the mean & stdev for the subset
	for ( int i = 0; i < binsize; i++ )
	{
		mean[ i ] = ComputeMean( row, 162000 / binsize );
		stdev[ i ] = ComputeStdev( row, 162000 / binsize, mean[ i ] );
		row += 81; //increment to the next set of 81
	}
	ComputedMean = ComputeMean( mean, binsize ); //cache this here so we don't have to compute it more than once
	ComputedStdev = ComputeStdev( mean, binsize, ComputedMean );

	//2c

	//Confidence intervals to use in checking number of times SIM is inside 90% CI range
	CILow = ComputedMean - ( 1.64*sqrt( 81 * pow( ComputedStdev, 2 ) ) ) / 9;
	CIHigh = ComputedMean + ( 1.64*sqrt( 81 * pow( ComputedStdev, 2 ) ) ) / 9;

	//Dumb code that computes the number of times the SIM falls in the 90% CI range
	for ( int i = 0; i < binsize; i++ )
	{
		if ( ( mean[ i ]> CILow ) && ( mean[ i ] < CIHigh ) )
		{
			NumTimesInsideCI++;
		}
	}

	//print the results
	printf( "Mean of SIMs: %-lf\n", ComputedMean );
	printf( "Stdev of SIMs: %-lf\n", ComputeStdev( mean, binsize, ComputedMean ) );
	printf( "Number of times SIMs inside CI: %d\n", NumTimesInsideCI );
} //end IntervalMeanAndDev

///<summary>
///Computes and outputs two CSV files (Histogram_<paramref name="rownumber"/>.csv, and HistogramPDF_ <paramref name="rownumber"/>.csv)
///</summary>
///<param name="rownumber">The row number to use in generating the filename</param>
void Histogram( double *row, int size, int rownumber )
{
	const int NUMBINS = 2000;
	const int LENGTH = (int)sqrt( size );
	const double PI = 4.0 * atan( 1.0 );
	char Filename[ 64 ];
	int Histogram[ NUMBINS ];
	double Bins[ NUMBINS ], PDF[ NUMBINS ], GaussPDF[ NUMBINS ];
	double Max, Min; //the max and min of the array
	double scale; //for the gauss stuff
	double computedMean, computedStdev;
	double ASD = 0, sum = 0;

	SearchForMaxMin( row, size, &Max, &Min );
	LoadHistogramFromVector( Histogram, LENGTH, row, size, Max, Min );
	ComputeHistogramBins( Bins, int( LENGTH ), Max, Min );
	sprintf( Filename, "Histogram%d.csv", rownumber ); //copy everything into a pointer to the filename, getting it ready for the WriteHistogram()
	WriteHistogram( Filename, Histogram, Bins, LENGTH ); //write out the raw data

	//Convert Histogram to PDF
	for ( int m = 0; m < LENGTH; m++ )
	{
		PDF[ m ] = (double)Histogram[ m ] / (double)size;
	}
	sprintf( Filename, "HistogramPDF%d.csv", rownumber );//copy everything into a pointer to the filename, getting it ready for the WriteHistogram()
	WritePDF( Filename, PDF, Bins, LENGTH ); //write out the PDF

	//Gauss PDF
	computedMean = ComputeMean( row, size ); //precompute stdev for efficiency
	computedStdev = ComputeStdev( row, size, computedMean );
	double Scale = ( Bins[ 1 ] - Bins[ 0 ] ) / ( sqrt( 2 * PI ) * computedStdev );
	for ( int m = 0; m < LENGTH; m++ )
	{
		GaussPDF[ m ] = Scale * exp( -( Bins[ m ] - computedMean ) *( Bins[ m ] - computedMean ) / ( 2.0*pow( computedStdev, 2 ) ) );
	}
	sprintf( Filename, "PDF%d.csv", rownumber );
	WritePDF( Filename, GaussPDF, Bins, LENGTH );

	//ASD
	for ( int i = 0; i < NUMBINS; i++ )
	{
		sum += fabs( GaussPDF[ i ] - PDF[ i ] );
	}
	printf( "ASD is %lf\n\n", sum );
} //end Histogram

///<Summary>
/// Function to write a histogram to a file.
/// Note the histogram is an array of integers, which contain the number of times
/// the data being histogrammed fell within a bin.
/// Also included with the counts, is a double precision number that is the
/// center value for each bin.
///</Summary>
void WriteHistogram( char *name, int *Histo, double *Bins, int bins )
{
	// Open the file.
	FILE *fout = fopen( name, "w" );
	// Check for valid file open.
	if ( fout )
	{
		// Loop through the data.
		while ( bins-- )
		{
			fprintf( fout, "%d,%18.16lg\n", // writeout and
				*Histo++, *Bins++ ); // move to next entries.
		} // End of loop through bins.
		fclose( fout );
	} // End of valid file open test.
} // End of WriteHistogram

/// <summary>Function to write a Probablity Distribution to a file.
/// This is similar to the histogram write about, except
/// instead of the integer counts, a double precision probablity is given.
/// </summary>
void WritePDF( char *name, double *Pdf, double *Bins, int bins )
{
	// Open the file.
	FILE *fout = fopen( name, "w" );
	// Check for valid file open.
	if ( fout )
	{
		// Loop through the data.
		while ( bins-- )
		{
			fprintf( fout, "%18.16lg,%18.16lg\n", // writeout and
				*Pdf, *Bins ); // move to next entries.
			Pdf++; Bins++;
		} // End of loop through bins.
		fclose( fout );
	} // End of valid file open test.
} // End of WriteHistogram

///<summary>The main method.</summary>
int main( )
{
	matrix InputMatrix; //the main matrix that we're reading in
	matrix row1, row2, row3, row4, row5; //separate matricies for each row
	matrix FaultTol1, FaultTol2, FaultTol3, FaultTol4;
	double CorrelCoeff1, CorrelCoeff2, CorrelCoeff3, CorrelCoeff4;

	//1a)

	InputMatrix = ReadBinaryMatrix( "StatsData.mtx" ); //read in the matrix from file
	//calculate the boundaries for the matricies so we don't have to calculate it multiple times
	const int NUMCOLS = InputMatrix.wide( );
	const int SAMPLESIZE = 81; //SampleSize for Confidence Interval

	//Allocate the space for the rows
	row1 = matrix( NUMCOLS, 1 );
	row2 = matrix( NUMCOLS, 1 );
	row3 = matrix( NUMCOLS, 1 );
	row4 = matrix( NUMCOLS, 1 );
	row5 = matrix( NUMCOLS, 1 );

	//matrix is created in this format (rows, cols). Don't forget row1 is InputMatrix row0
	for ( int i = 0; i < NUMCOLS; i++ )
	{
		row1( i ) = InputMatrix( 0, i );
		row2( i ) = InputMatrix( 1, i );
		row3( i ) = InputMatrix( 2, i );
		row4( i ) = InputMatrix( 3, i );
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

	printf( "Linear Regression (LR) Parameters & Correlation Coefficient (CC): \n" );
	printf( "LR: Row 1" );
	PrintMatrix( FaultTol1 );
	printf( "CC Row Row 1 & 5: %lg \n\n\n", CorrelCoeff1 );

	printf( "LR: Row 2" );
	PrintMatrix( FaultTol2 );
	printf( "CC Row Row 2 & 5: %lg \n\n\n", CorrelCoeff2 );

	printf( "LR: Row 3" );
	PrintMatrix( FaultTol3 );
	printf( "CC Row Row 3 & 5: %lg \n\n\n", CorrelCoeff3 );

	printf( "LR: Row 4" );
	PrintMatrix( FaultTol4 );
	printf( "CC Row Row 4 & 5: %lg \n\n\n", CorrelCoeff4 );

	//printf( "Correlation Coefficients\n" );

	//1b)

	//Row 3 and Row 4 have the highest CC, so perform a multivariable regression on them
	printf( "Multivariable Regression:\n" );
	MultiVariableRegression( row3, row4, row5, NUMCOLS );

	//2a)
	printf( "\n" );
	printf( "Row 1 stats:\n" );
	ConfidenceInterval( row1, SAMPLESIZE );
	printf( "Row 2 stats:\n" );
	ConfidenceInterval( row2, SAMPLESIZE );
	printf( "Row 3 stats:\n" );
	ConfidenceInterval( row3, SAMPLESIZE );
	printf( "Row 4 stats:\n" );
	ConfidenceInterval( row4, SAMPLESIZE );
	printf( "Row 5 stats:\n" );
	ConfidenceInterval( row5, SAMPLESIZE );

	//2b - d)
	printf( "SIMs:\n" );
	printf( "Row 1:\n" );
	IntervalMeanAndDev( row1.AsPointer( ) );
	Histogram( row1.AsPointer( ), NUMCOLS, 1 );
	printf( "Row 2:\n" );
	IntervalMeanAndDev( row2.AsPointer( ) );
	Histogram( row2.AsPointer( ), NUMCOLS, 2 );
	printf( "Row 3:\n" );
	IntervalMeanAndDev( row3.AsPointer( ) );
	Histogram( row3.AsPointer( ), NUMCOLS, 3 );
	printf( "Row 4:\n" );
	IntervalMeanAndDev( row4.AsPointer( ) );
	Histogram( row4.AsPointer( ), NUMCOLS, 4 );
	printf( "Row 5:\n" );
	IntervalMeanAndDev( row5.AsPointer( ) );
	Histogram( row5.AsPointer( ), NUMCOLS, 5 );

	//Generate histogram CSV files (both raw data and PDF)

	getchar( );

	return 0;
}