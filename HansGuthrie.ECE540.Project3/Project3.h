void MultiVariableRegression( matrix a, matrix b, matrix dependent, int size );

void ConfidenceInterval( matrix m );

void IntervalMeanAndDev( matrix *row );

void Histogram( double *row, int size, int rownumber );

void WriteHistogram( char *name, int *Histo, double *Bins, int bins );

void WritePDF( char *name, double *Pdf, double *Bins, int bins );

int main( );