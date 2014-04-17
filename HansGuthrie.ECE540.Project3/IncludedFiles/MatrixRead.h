#include "matrix.hpp"
#include <stdio.h>
//#include "ReadInCsvFile.h"

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
//Even though this a header file, we will code these functions here 
//since we want to use them in future projects

#ifndef MatrixRead_h
#define MatrixRead_h 0

////Reads in a matrix from a CSV file. 
//matrix ReadCsvMatrix( char *Name )
//{
//	matrix output;
//	int j, k, counter = 0;
//	int    Length = 0, Rows, Columns;
//	double *Data, temp;
//
//	Data = ReadInCSVFile( Name, &Rows, &Columns );
//	//matrix output( Rows, Columns );
//	output = matrix( Rows, Columns );
//
//	// If invalid read
//	if ( !Data )
//	{
//		// tell user and exit.
//		printf( "Error reading %s\n", Name );
//		return 0;
//	}
//
//	//loop through rows to transform the Data array into a matrix
//	for ( j = 0; j < Rows; j++ )
//	{
//		//loop through the columns
//		for ( k = 0; k < Columns; k++ )
//		{
//			//take this point (j, k) and copy the appropriate value of the Data array in to it
//			output( j, k ) = Data[ counter ];
//			counter++; //increment the counter
//		}
//	}
//	
//	return output;
//} // end ReadCsvMatrix

//Reads in a matrix file from a binary formatted matrix. 
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
