#include "matrix.hpp"
#include <stdio.h>
#include "MatrixRead.h"
#include "Stats_Support.h"
#include "MatrixOutputs.hpp"
#include <vector>

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

void part1( )
{
	matrix InputMatrix;
	InputMatrix = ReadBinaryMatrix( "StatsData.mtx" );
	int numrows = InputMatrix.high( );
	int numcols = InputMatrix.wide( );
	//matrix* rows[ 5 ];
	matrix* rows = new matrix[ 6 ];

	//allocate space in each 'row matrix'
	/*for ( int i = 0; i < numrows; i++ )
	{
		
	}*/
	//allocates space in each row matrix and copies the matrix from InputMatrix to the appropriate row matrix
	for ( int i = 1; i <= numrows; i++ )
	{
		rows[ i ] = matrix( 1, numcols ); //allocate space in each 'row matrix'
		
		//actually copy the values from InputMatrix to the 'row matrix'
		for ( int j = 0; j < numcols; j++ )
		{
			rows[ i ]( j ) = InputMatrix( i, j );
		}
	}

	//fault tolerance between Row1 & Row5
	PrintMatrix( FaultTolerantRegression( rows[ 1 ], rows[ 5 ] ) );

	//delete[] rows;
}

void main( )
{
	part1( ); 

}