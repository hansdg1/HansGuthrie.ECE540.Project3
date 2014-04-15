
#include "matrix.hpp"
#include <stdio.h>

using namespace std;

#define NAN_Local ((float)((1e300*1e300) * 0.0F))

////////////////////////////////////////////////////////////
// used in declaration of empty matrix.
matrix::matrix()
{ data = NULL;  // No data
  rows = 0;
  columns = 0;
}


////////////////////////////////////////////////////////////
// Used to create a matrix of a given size. 
matrix::matrix(int r, int c)
{
	// Allocate space for data and check for valid allocation.
	data = new double[r*c];
	if (data)
	{
		// if not worried about time, initialize data to zero
#ifdef FaultTolerant
		for (int i = 0; i<(r*c); i++)
			data[i] = 0.0;
#endif
		rows = r;  // load values into object.
		columns = c;

#if debug    
		printf("Creating matrix of %d by %d\n", r, c);
#endif

	}
#if debug
	else
	{
		printf("Error creating matrix with %d by %d\n", r, c);
		getchar();
	}
#endif
} // end of constructor

////////////////////////////////////////////////////////////
// Constructor for vector
matrix::matrix(int r)
{
	// Allocate space for data and check validity
	data = new double[r];
	if (data)
	{
		// if not concerned about run time initialize memory to zero.
#ifdef FaultTolerant
		for (int i = 0; i<r; i++)
			data[i] = 0.0;
#endif
		rows = r;  // load parameters into object
		columns = 1;  // since its a vector, only one column.

#if debug    
		printf("Creating vector of %d \n", r);
#endif

	}
#if debug
	else
	{
		printf("Error creating vector of %d \n", r );
		getchar();
	}
#endif
} // end of vector constructor.


////////////////////////////////////////////////////////////
// Constructor used when passing a matrix in to a function.
// also called a copy constructor.
matrix::matrix(matrix &m)
{ 
  // allocate space for data and check for validity.
  data = new double [m.high()*m.wide()];
  if( data )
  {
	  int k;
	  // set parameters in object.
	  rows = m.rows;
	  columns = m.columns;
	  // if not worried about run time, do a direct copy 
#ifdef FaultTolerant
	  int j;
	  for( k=0; k<m.high()*m.wide(); k++ )
	  	for( j=0; j<m.wide(); j++ )
			data[k*columns+j] = m(k,j); 
#else 
	  // copy data between the two arrays.
	  int length = rows*columns;
	  for (k = 0; k < length; k++)
		  data[k] = m.data[k];
#endif
#if debug    
		printf("Copy of matrix of %d by %d\n", m.high(), m.wide() );
#endif

  } // end of validity test
#if debug
  else
  {
	  printf("Error creating matrix with %d by %d\n", m.high(), m.wide());
	  getchar();
  }
#endif

} // end of copy constructor.


////////////////////////////////////////////////////////////
// Destructor used to remove an object.
matrix::~matrix()
{
  if( data )
  {

#if debug    
	  printf("Destroying matrix of %d by %d\n", rows, columns);
#endif

	delete [] data;  // free memory in data.
  }
} // end of destructor.


////////////////////////////////
// Overload of equal operator.

matrix matrix::operator =(matrix m)
{
 // if data holds data.
 if( data )     // first destroy the matrix pointed to by this instance.
 {delete [] data; // remove array
  rows = 0; // force rows and columns to zero.
  columns = 0;
 }
 if(  ! m.isValid() )
 {data = NULL;
  rows = 0;
  columns = 0;
 }
 else
 {
     // allocate space for data and check for validity
	 data = new double [m.high()*m.wide()];
	 if( data )
	 {
	  // move parameters over.
	  int k;
	  rows = m.rows;
	  columns = m.columns;
	  
	  // copy data over.
#ifdef FaultTolerant
	  int j;
	  for( k=0; k<m.high(); k++ )
	  	for( j=0; j<m.wide(); j++ )
			data[k*columns+j] = m(k,j); 
#else
	  int length = rows*columns;
	  for (k = 0; k < length; k++)
		  data[k] = m.data[k];
#endif

#if debug    
		printf("Copy of matrix of %d by %d\n", m.high(), m.wide() );
#endif
	 
	  } // end of valid memory
  	  else // if memory invalid 
	  {
		  data = NULL;
		  rows = 0;
		  columns = 0;
#if debug
		  printf("Error Coping matrix with %d by %d\n", m.high(), m.wide());
	          getchar();
#endif
	  } // end of invalid memory allocation.

  } // end of valid input matrix.

  return *this;
} // end of = overload,

////////////////////////////////
// Overload of add operator.
matrix matrix::operator +(matrix m2)
{int k;
 // check that matrices match in size.
 if(  m2.wide() != this->wide() 
   || m2.high() != this->high() )
 {matrix temp;
  return temp;  // return invalid matrix.
 }
 else // if operation valid.
 {
  // create matrix.
  matrix temp(m2.high(),m2.wide());
  // add elements together.
#ifdef FaultTolerant
  for (k = 0; k<m2.high(); k++)
      for (int j = 0; j<m2.wide(); j++)
	      temp(k, j) = this->data[k*columns + j] + m2(k, j);
#else
  for (k = 0; k<m2.high()*m2.wide(); k++)
      temp(k) = this->data[k] + m2.data[k];
#endif
  return temp;
 } // end of valid operation.
} // end of + overload.


////////////////////////////////
// Overload of subtract operator.
matrix matrix::operator -(matrix m2)
{int k;
 // Check that operation is valid.
 if(  m2.wide() != this->wide() 
   || m2.high() != this->high() )
 {matrix temp;
  return temp; // return invalid matrix
 }
 else  // Valid operation.
 {
  // create a matrix 
  matrix temp(m2.high(),m2.wide());
  // do difference of data arrays.
#ifdef FaultTolerant
  int j;
  for( k=0; k<m2.high(); k++)
     for( j=0; j<m2.wide(); j++)
	    temp(k,j) = this->data[k*columns+j] - m2(k,j);
#else
 int length = m2.rows*m2.columns;
 for (k = 0; k < length; k++)
	 temp.data[k] = data[k] - m2.data[k];
#endif
  return temp;
 } // end of valid operation
} // end of - overload.

////////////////////////////////
// Overload of matrix multiply.
matrix matrix::operator *(matrix m2)
{int k,j,m;
 // Check for valid operation.
 if( m2.high() != this->wide() )
 {matrix temp;
  return temp;  // if invalid return invalid matrix.
 }
 else // operation valid.
 {
  // create matrix matching size of result.
  matrix temp(this->high(),m2.wide());
  double tmp;
  // Perform operation.
#ifdef FaultTolerant
  for( k=0; k<this->high(); k++) // loop through rows.
     for( j=0; j<m2.wide(); j++) // loop through columns.
	 {
		tmp = 0;
		for( m=0; m < this->wide(); m++ )  // sum for element k,j.
			tmp += (*this)(k,m) * m2(m,j); 
		temp(k,j) = tmp;
	 }
#else
  // this is a pointer based solution for faster implementation.
  int length = temp.rows*temp.columns;
  double *t_ptr, *m2_ptr;
  for (k = 0; k < rows; k++)
  {
	  for (j = 0; j < m2.columns; j++)
	  {
	     m2_ptr = m2.data + j;
		 t_ptr = data + k*columns;
		 tmp = *t_ptr * (*m2_ptr);
		 for (m = 1; m < this->columns; m++)
		 {
			 t_ptr++;
			 m2_ptr += m2.columns;
			 tmp += *t_ptr * (*m2_ptr);
		 }
		 temp(k, j) = tmp;
	  }
  }
#endif
  return temp;
 } // end of valid operation check
} // end of * overload.

////////////////////////////////
// Overload of multiply by a scalar.
matrix matrix::operator *(double d)
{
	int k;
	matrix temp = *this;
#ifdef FaultTolerant
	int j;
	for( k=0; k < this->high(); k++)
	for( j=0; j < this->wide(); j++)
	{
		temp(k,j) = d*temp(k,j);
	}
#else
	int length = rows*columns;
	for (k = 0; k < length; k++)
		temp.data[k] = d*data[k];
#endif
	return temp;
} // end of scalar multiply////////////////////////////////

// Overload  Raise to a power operation
matrix matrix::operator^(int N)
{
	int k;
	matrix temp = *this;
	for (k = 1; k < N && temp.isValid(); k++)
		temp = temp * (*this);
	return temp;
} // end of Raise to a power operation

////////////////////////////////
// Overload of divide by a scalar.
matrix matrix::operator /(double d)
{
	double scale = 1.0 / d;
	return *this*scale;
} // end of divide by scalar.

double & matrix::operator()(int r, int c)
{
	// if you are not interested in run time, 
	// this will expand matrix if you are writing outside the matrix bounds.
#ifdef FaultTolerant
	if (r >= rows     // check for out of bounds
		|| c >= columns)
	{
		// access is outside of bounds.
		int R, C;
		R = (r >= rows) ? r + 1 : rows;  // look for new size
		C = (c >= columns) ? c + 1 : columns;
		double *local = new double[R*C];  // allocate space for new matrix.
		if (local)
		{
			// valid allocation.
#ifdef debug
			printf("Expanding matrix to %d by %d\n", R, C);
#endif
			// Loop through input matrix.
			for (int k = 0; k < R; k++)
			  for (int m = 0; m < C; m++)
			  {
				if (k < rows && m < columns)
					local[k*C + m] = this->data[k*columns + m];
				else
					local[k*C + m] = 0.0;
			  } // end of loop through input matrix.
			// Remove current array
			delete[](this->data);
			// Move data over to object 
			rows = R;
			columns = C;
			this->data = local;
		}// End of valid allocation test.
	} // end of test for out of bounds.
#endif
	// access array.
	return data[(r*columns + c)];
} // end of &operator()

/////////////////////////////////////////////////////
// () operator with only one index. For the case of a vector.
double & matrix::operator()(int index)
{
	// if access is outside of matrix.
	if (index < 0 || index >= rows*columns)
    {
		static double x = NAN_Local;
		return x;  // return not a number.
	}
    // return data 
	return data[index];

} // end of &operator()

////////////////////////////////////////////////////
// returns the transpose of the matrix
matrix matrix::transpose()
{

	matrix temp;
	if (isValid())
	{

		// Allocate space for matrix.
		temp = matrix(columns, rows);
		if (temp.isValid())
		{
			// loop through rows.
			for (int i = 0; i < rows; i++)
			{
				// loop through columns.
				for (int j = 0; j < columns; j++)
				{
					temp(j, i) = data[i*columns + j];
				}// end of loop through columns

			} // end of loop through rows.

		}// end of valid output check

	} // end of valid input check

	return temp;
} // end of transpose.

////////////////////////////////////////////////////
// returns the trace of the matrix
double matrix::trace()
{

	double temp = NAN_Local;
	if (isValid())
	{
		temp = 0.0;

		// Allocate space for matrix.
		for (int i = 0; i < rows && i < columns; i++)
			{
				temp += data[i*columns + i];

			} // end of loop through diagonal.

	} // end of valid input check

	return temp;
} // end of trace.

// returns the cofactor matrix that results from 
// eliminating a specified row and column
// Private function for internal use only.
matrix matrix::cofact(int row, int col) 
{
	matrix temp; // start with empty matrix.

	// Check for valid operation
	if( isValid()  // matrix is valid
		&& row < rows && rows >= 0 // input row is valid
		&& col < columns && col >= 0 )  // input col is valid.
	{  
	   temp = matrix(rows-1, columns-1);
       // copy the matrix, except for the removed row and column
	   for (int i = 0; i < this->rows - 1; i++)
	   {
		   for (int j = 0; j < this->columns - 1; j++)
		   {
			   if (i < row && j < col)
				   // before the deleted row and deleted column
				   temp(i, j) = this->data[i*columns + j];
			   else if (i < row && j >= col)
				   // before the deleted row, after the deleted column
				   temp(i, j) = this->data[i*columns + j + 1];
			   else if (i >= row && j < col)
				   // after the deleted row, before the deleted column
				   temp(i, j) = this->data[(i + 1)*columns + j];
			   else
				   // after the deleted row and column
				   temp(i, j) = this->data[(i + 1)*columns + j + 1];
		   }
	   }
	}
	return temp; // empty is returned if invalid operation.
}  // end of cofact()

// returns the determinant of the matrix, computed recursively
double matrix::det()
{
	if (this->rows == 0
		|| this->columns == 0)
	{	// verify matrix is valid
		return 0;
	}
	if (this->rows == 1
		&& this->columns == 1)
	{  // determinant of a single value is that value
		return data[0];
	}

	double D = 0, sign_flip = 1.0;
	// recursively compute the determinant
	for (int i = 0; i < rows; i++)
	{
		D += sign_flip * data[i*columns] * cofact(i, 0).det();
		sign_flip *= -1.0;
	}
	return D;
}

// returns the determinant of the matrix, computed using the SVD of the matrix. 
double matrix::fast_det()
{
	double d;
	matrix a = *this;
    // Check that matrix is square and valid.  
	if (a.wide() == a.high() && isValid() )
	{
		// find minimum of rows and columns
		int C = (rows > columns) ? rows : columns;

		// Allocate space for output, and check for validity
		matrix U, D(a.high(), a.wide()), V(a.wide(), a.wide());

		U = dsvd(a, D, V);

		// Check for valid allocations of matrix 
		if (U.isValid())
		{
			d = D(0,0);
			for (int k = 1; k < C; k++)
				d *= D(k, k);

		} // End of valid svd check.

	}// End of check for valid input matrix.

	return d; // return result, with Out being NULL if anything fails.

} // End of fast_det()

// returns the matrix inverse computed using the adjoint method
matrix matrix::inverse_det() {
	if(  this->rows == 0 
	  || this->columns == 0 ) {  // verify matrix is valid
		matrix temp;
		return temp;
	}

	matrix inv(this->rows, this->columns);
	double det_orig = det(), sign_flip = -1.0;

	for( int i = 0; i < this->rows; i++) {
		for( int j = 0; j < this->columns; j++) {
			inv(i, j) = sign_flip * this->cofact(i,j).det();
			inv(i, j) /= det_orig;
			sign_flip *= -1.0;
		}
	}
	
	// note transpose is performed outside the loop to make its use explicit
	inv = inv.transpose();	
	return inv;
}
// returns the matrix inverse computed using the LU inverse
matrix matrix::inverse_lu() 
{
	int i, j;
	// Check for valid matrix for inversion
	if(  this->rows == 0 || this->columns == 0 
	  || this->rows != this->columns) 
	{  // verify matrix is valid
		matrix temp;
		return temp;  // returns empty matrix on failure.
	}// end of valid matrix check.

    // Set up matrix and vectors for doing inverse.
	matrix inv(rows, columns), b, x(rows);
    matrix LU_A = *this;
    int *p = new int [rows];

	if(  inv.isValid() 
	  && b.isValid() 
	  && x.isValid()
	  && p )
	{
		// create lu factorization
	   LU( LU_A, p );

	   // move through columns
 	   for( i = 0; i < this->columns; i++) 
	   {
		 // create basis vector.
	  	 for( j = 0; j < this->rows; j++ )
			x(j) = 0.0;
		 x(i) = 1.0;

		 // solve for b = inv(A)*x;
		 b = usolve( LU_A, lsolve( LU_A, permutate( x, p ) ) );

		 // copy b to inverse
		 for( j = 0; j < this->rows; j++) 
			inv(j, i) = b(j);

	   }// end of loop through rows.

	}// end of validity check

	// delete array of indices.
	if (p) delete[] p;

	return inv;

}// end of inverse_det()
	
///////////////////////////////////////////////////////////////////
// This is a fault tolerant inversion functions which uses the SVD 
// to implement a matrix inverse, used in the psuedo inverse of a matrix.
matrix matrix::inverse( double Tolerance)
{
	matrix Out;  // Starts as a null matrix.
    double Max_SV;
	// Check that matrix is square and valid.  
	if (columns == rows && isValid())
	{
		// make a copy of a.
		Out = *this;

		// Allocate space for output, and check for validity
		matrix U, D(rows, columns), V(columns, columns);

		U = dsvd(*this, D, V);

		// Check for valid allocations of matrix 
		if (U.high() == rows)
		{
			int C = D.high();
			if (C > D.wide())
				C = D.wide();
			Max_SV = D(0);
			for (int k = 1; k < C; k++)
			    if( Max_SV < D(k,k) )
				    Max_SV = D(k,k);
			Max_SV = 1.0/Max_SV;
			
			for (int k = 0; k < C; k++)
				if( D(k,k)*Max_SV > Tolerance * eps )
				    D(k, k) = 1.0 / D(k, k);
				else
					D(k,k) = 0.0;

			Out = V*D*U.transpose();
		} // End of valid svd check.

	}// End of check for valid input matrix.

	return Out; // return result, with Out being NULL if anything fails.

} // End of inverse

///////////////////////////////////////////////////////////////////
// This is a fault tolerant inversion functions which uses the SVD 
// to implement a matrix inverse, used in the psuedo inverse of a matrix.
double matrix::condition()
{
	double d_max,d_min;  // Starts as a null matrix.

	// Check that matrix is square and valid.  
	if (columns == rows && isValid())
	{
		// Allocate space for svd output, and check for validity
		matrix U, D(rows, columns), V(columns, columns);

		U = dsvd(*this, D, V);

		// Check for valid allocations of matrix 
		if (U.isValid())
		{
			int C = D.high();
			if (C > D.wide())
				C = D.wide();
			d_max = D(0, 0);
			d_min = D(0, 0);
			for (int k = 1; k < C; k++)
			{
				if (d_max < D(k, k))
					d_max = D(k, k);
				if (d_min > D(k, k))
					d_min = D(k, k);
			}

		} // End of valid svd check.

	}// End of check for valid input matrix.

	return d_max/d_min; // return result.

} // End of inverse


// returns true if the matrix is valid, computations can be performed on it
bool matrix::isValid() {
	if (rows == 0 || columns == 0 || !data ) {  // verify matrix is valid
		return false;
	}

	return true;
} // end of isValid()

// returns true if the matrix is valid, computations can be performed on it
bool matrix::isInvalid() {
	if (rows == 0 || columns == 0 || !data) {  // verify matrix is valid
		return true;
	}

	return false;
} // end of isInvalid()

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
/ The LU function is to be used in conjunction with the 
/  usolve, lsolve and permutate to solve the system of equations.
/
/ To solve A x = b, would be as follows
  matrix A, LU_A, b, x;
  int p[rows];

  LU_A = A;
  LU( LU_A, p );
  x = usolve( LU_A, lsolve( LU_A, permutate( b, p ) ) );

*/
void LU(matrix &m, int *p) 
{
	// error check the matrix
	if( m.isInvalid()) 
	{
		return;
	}
	// error check the pointer
	if(p == NULL) 
	{
		return;
	}
	
	// initialize permutation vector
	for( int i = 0; i < m.high(); i++ )
		p[i] = i;

	// compute the LU factorization
	int pivot, temp2;
	double temp;
	for( int k = 0; k < m.wide() - 1; k++ )
	{

		// for pivoting determine if there is a 
		// value in the current column that has a 
		// larger magnitude than the current on diagonal value
		pivot = k;
		for( int i = k+1; i < m.high(); i++ ) 
		{
			if( abs(m(i,k)) > abs(m(pivot,k)) )
				pivot = i;
		}

		// pivot if necessary by exchaning rows
		if( pivot > k ) {
            for( int j = 0; j < m.wide(); j++ ) 
			{
				temp = m(k,j);
				m(k,j) = m(pivot, j);
				m(pivot, j) = temp;
			}  // end of loop through columns
			// update permutation vector
			temp2 = p[k];
			p[k] = p[pivot];
			p[pivot] = temp2;
		} // end of need for pivot check.

		// skip the current column if the diagonal is zero
		if( m(k, k) == 0 )
			continue;

		// compute the multipliers for elimination (L elements)
		for ( int i = k+1; i < m.high(); i++ ) 
		{
			m(i, k) = m(i,k)/m(k,k);	// store elements "in place"
		}
	
		// propagate the elimination through the rest of the matrix
		for( int j = k+1; j < m.wide(); j++ )
			for( int i = k+1; i < m.high(); i++)
				m(i,j) = m(i,j) - m(i,k)*m(k,j);  // store elements "in place"
	
	}
} // end of lu factorization

/*
/ Permute switches the rows of matrix m so that they are in the order
/ described by the order in p.  The number of rows in matrix m must 
/ match the number of elements in p.
/
/ inputs: m - matrix to permute
/		  p - new order of rows, same number of elements as rows in m 
*/
matrix permutate(matrix &m, int p[]) {
	// error checking
	if( p == NULL ) {
		matrix temp;
		return temp;
	}

	matrix temp(m.high(), m.wide());

	for( int k = 0; k < m.high(); k++ )
		for( int j = 0; j < m.wide(); j++ )
			temp(k,j) = m(p[k],j);
	
	return temp;
} // end of permutate

/*
/ lsovle solves the equation Ly=b using the L matrix from LU
/ factorization.
/
/ inputs: L - lower triangular L matrix from LU factorization, 
/				can be in the "in place" form with U in the 
/				upper diagonal
/		  b - solution matrix where each column is a solution
*/
matrix lsolve(matrix L, matrix b) {
	// error checking
	if( L.high() != b.high() ) {
		matrix temp;
		return temp;
	}

	for( int k = 0; k < L.high(); k++ )
		// loop over multiple solutions in b
		for( int j = 0; j < b.wide(); j++)	
			// update undetermined variables
			for( int i = k+1; i < L.high(); i++ )
				b(i, j) -= L(i,k)*b(k,j);
	
	return b;
} // end of lsolve.

/*
* usolve solves the equation Ux=y using the matrix U from LU
* factorization.
*
* inputs: U - Upper diagonal matrix from LU factorization.
*				can be in the "in place" form with L in the 
*				lower triangle.
*		  y - the solution matrix computed by lsolve()
*/
matrix usolve(matrix U, matrix y) {
	// error checking
	if( U.high() != y.high() ) {
		matrix temp;
		return temp;
	}

	// loop over columns, starting at far right
	for( int k = U.wide() - 1; k >= 0; k-- )
		// loop over multiple solution vectors in y
		for( int j = 0; j < y.wide(); j++ ) {  
			// set the next solution
			y(k,j) = y(k,j)/U(k,k);
			// update the undetermined variables based on 
			// the newly determined solution
			for(int i = k-1; i >= 0 ; i--)
				y(i,j) -= U(i,k)*y(k,j);
		}
	return y;
} // end of usolve

////////////////////////////////////////////////////////////////////////
// Support routines for SVD.
#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#define MAX(x,y) ((x)>(y)?(x):(y))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double PYTHAG(double a, double b)
{
    double at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
} //end of PYTHAG
// end of support routines.

/*-----------------------------------------------------------------------------
 The following routine computes the SVD of the matrix In and returns 
  U as the return value and is the left vectors,
  D as a diagonal matrix 
  V as the Right Vectors
  Note: D and V are assumed to be of the size needed to hold the data.
  if In is m by n, then U is an 
*/
// Please note that this is old translated code and is not written well,
// but it is working so don't touch it.
matrix dsvd(matrix In, matrix &D, matrix &v )
{
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;

	int m = In.high(),
		n = In.wide();

    matrix a;

    // Restrict high > wide. (If not the case (simply use transpose of matrix )
    if (m < n)
    {
        return a; // Return empty matrix as indicator of fault.
    }

    a = In;
    
	// Allocate work space 
	matrix w(m+n);  
    matrix rv1(m+n);

    /* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1(i) = scale * g;
        g = s = scale = 0.0;
        if (i < m) 
        {
			// Test size of column vector(k:m,i)
            for (k = i; k < m; k++)
                scale += fabs( a(k,i));
            // if not too small
			if (scale) 
            {
				// Rescale vector
                for (k = i; k < m; k++) 
                {
                    a(k,i) =  ( a(k,i)/scale);
                    s += ( a(k,i) *  a(k,i));
                }

				// Compute transform 
                f =  a(i,i);
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a(i,i) =  (f - g);
                if (i != n - 1)
                {
					// Apply transform
                    for (j = l; j < n; j++) 
                    {
                        for (s = 0.0, k = i; k < m; k++) 
                            s += ( a(k,i) *  a(k,j));
                        f = s / h;
                        for (k = i; k < m; k++) 
                            a(k,j) +=  (f *  a(k,i));
                    }
                }
                for (k = i; k < m; k++) 
                    a(k,i) =  ( a(k,i)*scale);
            }
        }
        w(i) =  (scale * g);
    
        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1)
        {
            for (k = l; k < n; k++) 
                scale += fabs( a(i,k));
            if (scale) 
            {
                for (k = l; k < n; k++) 
                {
                    a(i,k) = a(i,k)/scale;
                    s += ( a(i,k) * a(i,k) );
                }
                f = a(i,l);
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a(i,l) = (f - g);
                for (k = l; k < n; k++)
                    rv1(k) = a(i,k) / h;
                if (i != m - 1) 
                {
                    for (j = l; j < m; j++) 
                    {
                        for (s = 0.0, k = l; k < n; k++) 
                            s += ( a(j,k) * a(i,k) );
                        for (k = l; k < n; k++)
                            a(j,k) += (s * rv1(k));
                    }
                }
                for (k = l; k < n; k++) 
                    a(i,k) = ( a(i,k) * scale );
            }
        }
        anorm = MAX(anorm, (fabs(w(i)) + fabs(rv1(i))));
    }

    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        if (i < n - 1) 
        {
            if (g) 
            {
                for (j = l; j < n; j++)
                    v(j,i) = ( ( a(i,j) / a(i,l) ) / g );
                    /* double division to avoid underflow */
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < n; k++) 
                        s += ( a(i,k) * v(k,j) );
                    for (k = l; k < n; k++)
                        v(k,j) += (s * v(k,i));
                }
            }
            for (j = l; j < n; j++) 
                v(i,j) = v(j,i) = 0.0;
        }
        v(i,i) = 1.0;
        g = rv1(i);
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        l = i + 1;
        g = w(i);
        if (i < n - 1) 
            for (j = l; j < n; j++)
                a(i,j) = 0.0;
        if (g) 
        {
            g = 1.0 / g;
            if (i != n - 1) 
            {
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += ( a(k,i) * a(k,j) );
                    f = ( s / a(i,i) ) * g;
                    for (k = i; k < m; k++) 
                        a(k,j) += (f * a(k,i) );
                }
            }
            for (j = i; j < m; j++) 
                a(j,i) = ( a(j,i) * g );
        }
        else 
        {
            for (j = i; j < m; j++) 
                a(j,i) = 0.0;
        }
        ++a(i,i);
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--) 
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++) 
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1(l)) + anorm == anorm) 
                {
                    flag = 0;
                    break;
                }
                if (fabs( w(nm) ) + anorm == anorm) 
                    break;
            }
            if (flag) 
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1(i);
                    if (fabs(f) + anorm != anorm) 
                    {
                        g = w(i);
                        h = PYTHAG(f, g);
                        w(i) = h; 
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = a(j,nm);
                            z = a(j,i);
                            a(j,nm) = (y * c + z * s);
                            a(j,i) = (z * c - y * s);
                        }
                    }
                }
            }
            z = w(k);
            if (l == k) 
            {                  /* convergence */
                if (z < 0.0) 
                {   /* make singular value nonnegative */
                    w(k) =  (-z);
                    for (j = 0; j < n; j++)
                        v(j,k) = (-v(j,k));
                }
                break;
            }
            if (its >= 30) 
			{matrix tmp;
                return tmp;
            }
    
            /* shift from bottom 2 x 2 minor */
            x = w(l);
            nm = k - 1;
            y = w(nm);
            g = rv1(nm);
            h = rv1(k);
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
          
            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++) 
            {
                i = j + 1;
                g = rv1(i);
                y = w(i);
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1(j) = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++)
                {
                    x = v(jj,j);
                    z = v(jj,i);
                    v(jj,j) = (x * c + z * s);
                    v(jj,i) = (z * c - x * s);
                }
                z = PYTHAG(f, h);
                w(j) = z;
                if (z) 
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++) 
                {
                    y = a(jj,j);
                    z = a(jj,i);
                    a(jj,j) = (y * c + z * s);
                    a(jj,i) = (z * c - y * s);
                }
            }
            rv1(l) = 0.0;
            rv1(k) = f;
            w(k)   = x;
        }
    }

	// move singular values over to D matrix.
    for( k = 0; k < m; k++)
		 for( i = 0; i < n; i++ )
			 if( k == i )
				 D(k,i) = w( k );
			 else
				 D(k,i) = 0.0;

	return a;

}// end of dsvd

///////////////////////////////////////////////////////////////////
// Program that computes the pseudo inverse of a matrix.
matrix MatrixPseudoInverse( matrix a )
{
    matrix mcp, Imcp, Out;

	// Compute the (At * A) matrix.
    mcp = a.transpose()*a;
    
    if( mcp.high() > 0 ) // Check for valid
    {
		// Compute the inverse of(At * A)  
		Imcp = mcp.inverse();

        if( Imcp.high() > 0 ) // Check for valid 
        {
  		    // finally  inv(At * A) * At  or Pseudo Inverse of A.
            Out = Imcp*a.transpose();

		} // End of valid inverse mcp check

    } // End of valid mcp check

    return Out; // return result, with Out being NULL if anything fails.

} // End of MatrixPseudoInverse


// Quick Sort of a MtrxArray
// The following three routines are supporting functions
// for the main function named MatrixQuickSortColumn.
void MatrixRowSwitch(matrix *a, int pvt, int bot)
{
	register int k;
	double HoldValue;

	for (k = 0; k < a->wide(); k++)
	{
		HoldValue = (*a)(pvt, k);
		(*a)(pvt, k) = (*a)(bot, k);
		(*a)(bot, k) = HoldValue;
	}
}

int MatrixQSPartionCol(matrix *a, int Column, int top, int bot, bool Ascending)
{
	register int store, pvt;
	register double PvtValue;

	if (top >= bot)
		return bot;

	pvt = (top + bot) / 2;
	PvtValue = (*a)(pvt, Column);

	MatrixRowSwitch( a, pvt, bot);

	for (store = top;
		top < bot;
		top++ )
	{
		if (  ( Ascending && ((*a)(top, Column) < PvtValue))
		   || (!Ascending && ((*a)(top, Column) > PvtValue)))
		{
			if (top != store)
				MatrixRowSwitch( a, top, store);
			store++;
		}

	}
	if (store != bot)
	{
		MatrixRowSwitch( a, store, bot);
	}

	return store;

}

void MatrixQSRecurseCol(matrix *a, int Column, int top, int bot, bool Ascending)
{
	register int pvt;

	if (top < bot)
	{
		pvt = MatrixQSPartionCol(a, Column, top, bot, Ascending);
		MatrixQSRecurseCol(a, Column, top, pvt - 1, Ascending);
		MatrixQSRecurseCol(a, Column, pvt + 1, bot, Ascending);
	}
}
//End of support routines.

///////////////////////////////////////////////////////////////////
// This function will sort the columns in the matrix such that column
// "Column" is in acsending order. Ascending will determine if the data 
// is in Ascending or Descending.
//
// The inputs are a matrix "a", and the index of the column "Column"
// that is to be sorted.  Since the data in the MtrxArray
// is accessed by reference, the data that is in a is what is
// sorted, so no return is needed.  However, it should be noted 
// that the data in "a" will be permanently changed.
void MatrixQuickSortColumn(matrix *a, int Column, bool Ascending )
{
	MatrixQSRecurseCol(a, Column, 0, a->high() - 1, Ascending);
} // end of MatrixQuickSortColumn

// The following three functions are support routines for the
// main sort row function named MatrixQuickSortRow
void MatrixColumnSwitch(matrix *a, int left, int store)
{
	double HoldValue;
	for (int k = 0; k < a->high(); k++)
	{
		HoldValue = (*a)(k, left);
		(*a)(k, left) = (*a)(k, store);
		(*a)(k, store) = HoldValue;
	}
}
int MatrixQSPartionRow(matrix *a, int Row, int left, int right, bool Ascending)
{
	register int store, pvt;
	register double PvtValue;

	if (left >= right)
		return right;


	pvt = (left + right) / 2;
	PvtValue = (*a)(Row, pvt);

	MatrixColumnSwitch( a, pvt, right);

	for (store = left;
		left < right;
		left++)
	{
		if (   ( Ascending && ((*a)(Row, left) < PvtValue))
			|| (!Ascending && ((*a)(Row, left) > PvtValue)))
		{
			if (left != store)
				MatrixColumnSwitch( a, left, store);
			store++;
		}

	}
	if (store != right)
	{
		MatrixColumnSwitch( a, store, right);
	}

	return store;

}

void MatrixQSRecurseRow(matrix *a, int Row, int left, int right, bool Ascending)
{
	register int pvt;

	if (left < right)
	{
		pvt = MatrixQSPartionRow(a, Row, left, right, Ascending);
		MatrixQSRecurseRow(a, Row, left, pvt - 1, Ascending);
		MatrixQSRecurseRow(a, Row, pvt + 1, right, Ascending);
	}
}
// End of support routines.

///////////////////////////////////////////////////////////////////
// This function will sort the rows in the matrix such that row
// "Row" is in order. Ascending will determine if the data is in
// Ascending or Descending.
//
// The inputs are a matrix "a", and the index of the row "SRow"
// that is to be sorted.  Since the data in the MtrxArray
// is accessed by reference, the data that is in a is what is
// sorted, so no return is needed.  However, it should be noted 
// that the data in "a" will be permanently changed.

void MatrixQuickSortRow(matrix *a, int Row, bool Ascending )
{
	MatrixQSRecurseRow(a, Row, 0, a->wide() - 1, Ascending);
}// End of MatrixQuickSortRow

bool MatrixReverseRows(matrix *a)
{
	int k, n;

	if (a->isValid())
	{
		// loop through rows.
		for (k = 0, n = a->high() - 1; k < n; k++, n--)
		{
			MatrixRowSwitch(a, k, n);
		} // end of loop through rows.
		return true;
	}// end of validity check
	return false;
}// end of MatrixReverseRows

bool MatrixReverseColumns(matrix *a)
{
	int k, n;

	if (a->isValid())
	{
		// loop through rows.
		for (k = 0, n = a->wide() - 1; k < n; k++, n--)
		{
			MatrixColumnSwitch(a, k, n);
		} // end of loop through rows.
		return true;
	}// end of validity check
	return false;
}// end of MatrixReverseColumns

// This function will create a matrix of random values of size r,c
// with optional parameters of standard deviation and mean.
// Mean and Std default to 0.0 and 1.0
matrix RandomMatrix(int r, int c, double STD, double Mean)
{
	// Create new matrix
	matrix A(r, c);
	// Check that it is valid
	if (A.isValid())
	{
		// loop through rows
		for (int k = 0; k < r; k++)
		{
			// loop through columns
			for (int m = 0; m < c; m++)
			{
				// generate randome number.
				A(k, m) = STD * GaussianRandomNumbers() + Mean;
			} // end of loop through columns
		}// end of loop through rows.
	} // end of validity check.
	return A;
} // end of Random Matrix 

// This function will create a vector of random values of size r
// with optional parameters of standard deviation and mean.
// Mean and Std default to 0.0 and 1.0
matrix RandomVector(int N, double STD, double Mean)
{
	// Create new matrix
	matrix A(N, 1);
	// Check that it is valid
	if (A.isValid())
	{
		// loop through rows
		for (int k = 0; k <N; k++)
		{
			// generate randome number.
			A(k) = STD * GaussianRandomNumbers() + Mean;
		}// end of loop through rows.
	} // end of validity check.
	return A;
} // end of Random vector

///////////////////////////////////////////////////////////////////
// Support Routine: Accumulates the sum of the cross product of two sub-vectors.  Meant to support a matrix multiply.
static double MtrxAccumProd(double *a, int ac, int aStride, double *b, int bStride)
{
	register double res;

	if (!ac)
		return 0.0;

	res = (*a) * (*b); // Initialize res with first element.

	while (--ac)
	{
		a += aStride;
		b += bStride;
		res += (*a) * (*b);
	} // End of loop through data.

	return res; // Return of results.

} // End of MtrxAccumProd

///////////////////////////////////////////////////////////////////
// Employs the Conjugant Gradient method to solve A x = b.
matrix MatrixConjugateGradient(matrix &A, matrix b)
{
	register int i, k, rows, columns;

	register double scale, Rdot, RdotOld;

	matrix x, r, d, tmp;
	register double *Aptr;

	x = RandomVector( A.high() );
	r =   matrix( A.high() );
	d =   matrix( A.high() );
	tmp = matrix( A.high() );

	if (  x.isValid() && r.isValid() 
	   && d.isValid() && tmp.isValid() )
	{
		rows = A.high();
		columns = A.wide();
		Aptr = A.AsPointer();
		for (i = 0; i < rows; i++)
		{
			r(i) = b(i)  - MtrxAccumProd(Aptr, columns, 1, x.AsPointer(), 1);
			d(i) = r(i);
			Aptr += columns;
		}

		Rdot = MtrxAccumProd(r.AsPointer(), rows, 1, r.AsPointer(), 1);

		for (k = 0; k < rows; k++)
		{
			tmp = A*d;
			scale = Rdot / MtrxAccumProd(d.AsPointer(), rows, 1, tmp.AsPointer(), 1);

			x = x + (d*scale);   //  MtrxLineAdd(x->data, rows, 1, d, 1, scale);
			r = r - (tmp*scale); //  MtrxLineAdd(r, rows, 1, tmp, 1, -scale);

			RdotOld = Rdot;
			Rdot = MtrxAccumProd(r.AsPointer(), rows, 1, r.AsPointer(), 1);

			d = d * ( Rdot / RdotOld );
			d = d + r; // MtrxLineAdd(d, rows, 1, r, 1, 1.0);

		}

	}

	return x;

}

///////////////////////////////////////////////////////////////////
// The following routine will do a simplex reduction on the matrix A.
int MtrxSimplexReduction( matrix *A)
{
	register int k, m, offset, flag;

	register int    PvtColumn, PvtRow;
	register double PvtColumnValue, PvtRowValue;

	register double scale;

	flag = 1;

	{
		// Search bottom Row for Pvt Column.
		offset = (A->high() - 1) * A->wide();

		// Note the first column is the limits so it is skipped. 
		PvtColumnValue = (*A)(A->high()-1,1);
		PvtColumn = 1;

		for (k = 2; k < A->wide(); k++)
		{
			if (PvtColumnValue > (*A)(A->high()-1,k) )
			{
				PvtColumnValue = (*A)(A->high()-1, k);
				PvtColumn = k;
			}
		}

		// Check if True pivot column was found.
		if (PvtColumnValue < 0.0)
		{

			// Search for Pivot Row, by testing ratio.
			for (k = 0; k < A->high() - 1; k++)
			{
				scale = (*A)(k,PvtColumn);
				if (fabs(scale) > eps)
				{
					scale = (*A)(k, 0) / scale;
					if (flag && scale > eps)
					{
						PvtRowValue = scale;
						PvtRow = k;
						flag = 0;
					}
					else if (scale > eps && scale < PvtRowValue)
					{
						PvtRowValue = scale;
						PvtRow = k;
					}

				}
			}

			// If a valid pivot row was found
			if (flag == 0)
			{

				// Normalize Pvt Row.
				for (k = 0; k<A->wide(); k++)
				{
					if (k != PvtColumn)
						(*A)(PvtRow, k) /= (*A)(PvtRow, PvtColumn);

				}

				(*A)(PvtRow, PvtColumn) = 1.0;

				// Reduce other rows using pivot row.
				for (k = 0; k<A->high(); k++)
				{
					if (k != PvtRow)
					{
						scale = (*A)(k,PvtColumn)
							/ (*A)(PvtRow,PvtColumn);

						for (m = 0; m<A->wide(); m++)
						{
							(*A)(k, m) = (*A)(k,m) - scale *(*A)(PvtRow, m);

						}
					}
				}
			}

		}

	}
	
	return (flag == 0);

} // end of MtrxSimplexReduction

// Create an N by N identity matrix.
matrix eye(int N)
{
	return eye(N, N);

} // Create identity Matrix

// Create an N by M identity matrix.
matrix eye(int N, int M)
{
	matrix Out(N, M);
	if (Out.isValid())
	{
		for (int k = 0; k < N; k++)
		{
			for (int m = 0; m < M; m++)
			{
				Out(k, m) = (k == m) ? 1.0 : 0.0;
			}
		}
	}

	return Out;

} // Create identity Matrix

