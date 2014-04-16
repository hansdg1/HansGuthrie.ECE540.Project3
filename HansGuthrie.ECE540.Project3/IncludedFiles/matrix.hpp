#ifndef MATRIX_H
#define MATRIX_H

// This program is intended to establish a matrix object.
// The object will have the prime operators and 

#include <math.h>
#include "RandomNumbers.h"

const double eps = 2.22e-16;

#ifndef NULL
#define NULL 0
#endif

//#define debug 1

// This class defines a basic matrix, upon which we will build 
// a rudimentry system. 

class matrix
{
private:
	    double *data;    // This is the data in the matrix,
	    int rows, columns;   // and these are the size of the matrix.
		matrix cofact(int r, int c);
public:
	    // basic constructor
	    matrix(int r, int c); // Creates a matrix of size N by M
	    matrix(int N);        // Creates a vector of length N (actually an N,1 matrix.

		// empty constructor
		matrix( );

		// Copy Constructor
        matrix( matrix &m);

		~matrix();
		
		// Basic data access.
		int wide() {return columns;}
		int high() {return rows;}
		double *AsPointer() { return data; }

		// operator overloads.
		matrix operator=(matrix m);
		matrix operator+(matrix m);
		matrix operator-(matrix m);
		matrix operator*(matrix m);
		matrix operator*(double a);
		
		matrix operator/(double a);
		matrix operator^(int n);

        // Access to data  
		double &operator()(int r, int c); // Allows user to read or write to location r,c
		double &operator()(int index);    // Allows user to read or write to vector location index.

		// Check if the matrix is or is not Valid of Invalid
		bool isValid();
		bool isInvalid();
		
		// Computes determinate via
		double det();      // Adjoint 
		double fast_det(); // SVD

		// Computes the condition of the matrix. 
		double condition();  // the closer to one the more invertible.

		matrix inverse_det(); // Inverse by determinate approach
		matrix inverse_lu();  // Inverse by LU factorization
		matrix inverse(double Tolerance = 1.0 );     // Inverse by SVD 
		matrix transpose();   // Transpose of a matrix.
		double trace(); // Trace of matrix.

};

////////////////////////////////////////////////////////////////////////
/*-----------------------------------------------------------------------------
 The following routine computes the SVD of the matrix In and returns 
  U as the return value and is the left vectors,
  D as a diagonal matrix 
  V as the Right Vectors
  Note: D and V are assumed to be of the size needed to hold the data.
  if In is m by n, then U is an m by m, D is m by n and V is n by n. 
*/
matrix dsvd(matrix In, matrix &D, matrix &v );

/* 
/ LU computes the LU factorization of the matrix m "in place".
/ The elements of m below the diagonal are replaced with the L
/ matrix elements and the upper diagonal elements of m are 
/ replaced with the U elements.  Pivoting is used to compute
/ the factorization and the permutation of the rows are stored
/ in the array pointed to by p.
/
/ inputs: m - reference to the matrix on which LU factorization
/			  will be performed.
/		  p - pointer to array of integers that will hold the 
/			  permutation vector.
/
/ references: 
*/
void LU(matrix &m, int *p); 

/*
/ Permutate switches the rows of matrix m so that they are in the order
/ described by the order in p.  The number of rows in matrix m must 
/ match the number of elements in p.
/
/ inputs: m - matrix to permute
/		  p - new order of rows, same number of elements as rows in m 
*/
matrix permutate(matrix &m, int p[]);

/*
/ lsovle solves the equation Ly=b using the L matrix from LU
/ factorization.
/
/ inputs: L - lower triangular L matrix from LU factorization, 
/				can be in the "in place" form with U in the 
/				upper diagonal
/		  b - solution matrix where each column is a solution
*/
matrix lsolve(matrix L, matrix b); 

/*
* usolve solves the equation Ux=y using the matrix U from LU
* factorization.
*
* inputs: U - Upper diagonal matrix from LU factorization.
*				can be in the "in place" form with L in the 
*				lower triangle.
*		  y - the solution matrix computed by lsolve()
*/
matrix usolve(matrix U, matrix y);

//////////////////////////////////////
// This function will compute the least square pseudo inverse
// of a matrix.
matrix MatrixPseudoInverse( matrix a );
matrix MatrixConjugateGradient(matrix &A, matrix b);

// Functions that will sort a matrix based on a column or row.
// Always in ascending order.
void MatrixQuickSortColumn(matrix *a, int Column, bool Ascending = true );
void MatrixQuickSortRow(matrix *a, int Row, bool Ascending = true );

bool MatrixReverseRows(matrix *a);
bool MatrixReverseColumns(matrix *a);

// Create matrices of random numbers.
matrix RandomMatrix(int, int, double Std = 1.0, double Mean = 0.0 );
matrix RandomVector(int N, double Std = 1.0, double Mean = 0.0);
// Create identity matrices
matrix eye(int N);
matrix eye(int N, int M);

// Perform Simplex Step 
int MtrxSimplexReduction(matrix *A);

#endif
