#ifndef RandomNumbers_H_

#define RandomNumbers_H_

void SetUniformSeed(unsigned int seed);

// This function will generate Gaussian random numbers,
// with a mean of zero and a variance of 1.
double GaussianRandomNumbers();

// Function to search through vector and find its maximum and minimum.
void SearchForMaxMin(double *vector, int length,
                     double *Maximum, double *Minimum);

// This function will compute the mean of the data in vector.
double ComputeMean(double *vector, int length);

// This function will compute the Standard deviation of the 
// data in vector.
double ComputeStdev(double *vector, int length, double Mean);

// This function will load the array Histogram, of size "bins"
// with a count of the number of samples in random that fall within
// the range defined in steps of ( ( Maximum - Minimum ) / bins ). 
void LoadHistogramFromVector(int *Histogram, int bins,  // Histogram
                             double *random, int length, // Vector of data
                             double Maximum, double Minimum); // Range.

#endif  // RandomNumber_H_