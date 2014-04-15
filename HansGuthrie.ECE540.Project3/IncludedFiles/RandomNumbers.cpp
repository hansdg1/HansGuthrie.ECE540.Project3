#ifndef RandomNumbers_H_

#define RandomNumbers_H_

#include <math.h>
#include <stdlib.h>

// Define of special constant.
#define TwoPi 6.283185307179586476

// Workspace used to generate Gaussian random variables. 
int SecondRandomNumber = 0;
double Uniform1, Uniform2;

void SetUniformSeed( unsigned int seed )
{
    srand( seed );
}

// This function will generate Gaussian random numbers,
// with a mean of zero and a variance of 1.
double GaussianRandomNumbers( )
{
    // if on second pass, simply use previously
    // found Uniform numbers.
    if( SecondRandomNumber )
    {
        SecondRandomNumber = 0;  // flag to restart next pass.

        return Uniform1*sin(Uniform2); // Return second random number.
    }
    else  // on first pass, generate two uniform numbers.
    {
        SecondRandomNumber = 1;
        Uniform1 = (double)( ( rand() % 32768 ) + 1 ) / 32768.0;
        Uniform2 = (double)( ( rand() % 32768 ) + 1 ) / 32768.0;

        // Perform first part of transformation.
        Uniform1 = sqrt( -2.0 * log( Uniform1 ) );
        Uniform2 = TwoPi * Uniform2;
        
    } // End of first/second pass if statement.

    return Uniform1*cos(Uniform2); // Return first random number.

} // End of GaussianRandomNumbers

// Function to search through vector and find its maximum and minimum.
void SearchForMaxMin( double *vector,  int length, 
                      double *Maximum, double *Minimum )
{
    *Maximum = *vector; // Initialize data to the first 
    *Minimum = *vector; // entry in the vector.

    // Loop through vector
    while(--length)
    {
        vector++; 
        if( *Maximum < *vector )  // if Maximum < vector data
            *Maximum = *vector;   //    replace Maximum
        if( *Minimum > *vector )  // if Minimum > vector data
            *Minimum = *vector;   //    replace Minimum
    } // End of loop through vector;

} // End of SearchForMaxMin

// This function will compute the mean of the data in vector.
double ComputeMean( double *vector, int length )
{
    double Sum, Scale;

    Scale = 1.0 / length; // compute scale for mean calculations.

    Sum = *vector; // Initialize to first entry.

    // Loop through vector.
    while( --length )
    {
        vector++;
        Sum += *vector;
    } // End of loop through data.

    return Scale * Sum; // Return mean

} // End of ComputeMean

// This function will compute the Standard deviation of the 
// data in vector.
double ComputeStdev( double *vector, int length, double Mean )
{
    double Sum, Scale;

    Scale = 1.0 / (length-1); // compute scale for std calculations.

    Sum = (*vector-Mean)*(*vector-Mean); // Initialize to first entry.

    // Loop through vector.
    while( --length )
    {
        vector++;
        Sum += (*vector-Mean)*(*vector-Mean);
    } // End of loop through data.

    return sqrt( Scale * Sum ); // Return std
} // End of ComputeStdev

// This function will load the array Histogram, of size "bins"
// with a count of the number of samples in random that fall within
// the range defined in steps of ( ( Maximum - Minimum ) / bins ). 
void LoadHistogramFromVector( int *Histogram,  int bins,  // Histogram
                              double *random,  int length, // Vector of data
                              double Maximum, double Minimum ) // Range.
{
    int index; // Workspace.

    // Initialize Histogram to zero.
    for( index = 0; index < bins; index++ )
    {
        Histogram[index] = 0;
    } // End of Histogram initialization.

    // Loop through vector of data.
    while( --length )
    {
        // Compute index of where this data should be in histogram.
        index = (int) ( bins * ( *random - Minimum ) 
                             / ( Maximum - Minimum ) );
        // Check to make sure index is in range.
        if( index >= bins )
            index = bins-1;
        if( index < 0 )
            index = 0;
        // Increment that entry in histogram.
        Histogram[index]++;

        random++; // Move to next entry in vector.

    } // End of loop through vector.
} // End of LoadHistogramFromVector

#endif  // RandomNumber_H_