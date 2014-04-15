// Header for Matrix output functions, one to the display
// and the another to a file.
#ifndef MatrixOutputs_H_
#define MatrixOutputs_H_ 

#include "Matrix.hpp"
#include "stdio.h"

// Simple Print Out of a Matrix.
void PrintMatrix( matrix IA, int cols = -1 )
{
    int i = 0, i2, k, m;
    
    if( cols == -1 ) // If cols is not provided.
        cols = IA.wide(); // set cols to default value.
    
    // Move through matrix, “cols” at at time.
    while( i < IA.wide() )
    {
       i2 = i+cols; // Starting at i and going to i2-1.
    
       if( i2 > IA.wide() ) // being sure that i2-1 less than
           i2 = IA.wide(); // number of columns in matrix.
    
       printf( "\n Columns %d to %d \n", i+1, i2 );
    
       // Loop through rows.
       for( k = 0; k < IA.high(); k++ )
       { 
           printf( "| ");
     
           for( m = i; m < i2; m++ )
               printf( "%12.6f, ", IA(k,m) );
           if( i2 != IA.wide() )
               printf( " ... \n" );
           else
               printf( " | \n" );
       } // End of rows loop.
    
       i = i2;
    
    } // Loop through columns.
    
    printf( "\n" );

} // End of PrintMtrxArray

// Write Matrix to a file.
void WriteMatrixToAFile( matrix A, char *name )
{
    int i, k;

    FILE *fout;

    // Open file and test for valid open.
    fopen_s( &fout, name, "w" );
    if( fout )  
    {
        // Move through the rows in the matrix.
        for( i = 0; i < A.high(); i++ )
        {
            // Loop through columns, but not the last point.
            for( k = 0; k < A.wide()-1; k++ )
            {
                // Write out this element and move to next element.
                fprintf( fout, "%18.16lg, ", A( i, k ) );
            }

            // Write last element in  with a new line instead of a comma.
            fprintf( fout, "%18.16lg\n", A( i, k ) );

        } // End of loop through rows.

        fclose( fout ); // Close file.

    } // End of Valid file open check.

} // End of WriteMtrxArray


#endif
