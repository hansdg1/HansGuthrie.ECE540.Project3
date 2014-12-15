% Remove variable and plots.
clear
close all

% Read in the data
load('HistogramPDF1.csv');
load('PDF1.csv');
load('HistogramPDF2.csv');
load('PDF2.csv');
load('HistogramPDF3.csv');
load('PDF3.csv');
load('HistogramPDF4.csv');
load('PDF4.csv');
load('HistogramPDF5.csv');
load('PDF5.csv');


%% 
% for row 1
% Plot the histogram estimate of the Histogram
figure(1);
bar( HistogramPDF1(:,2), HistogramPDF1(:,1) );
hold on

% Plot the pdf's using the Gaussian distribution of 90% CI
plot( PDF1(:,2), PDF1(:,1), 'r' );

title( 'Histograms and PDF''s ' );
xlabel( 'Random Variable (Unitless)' );
ylabel( 'Probability Density' );
grid
legend( 'Histogram', 'Mean PDF');
print -dpng Figure1.png
%% 
% for row 2
% Plot the histogram estimate of the Histogram
figure(2);
bar( HistogramPDF2(:,2), HistogramPDF2(:,1) );
hold on

% Plot the pdf's using the Gaussian distribution of 90% CI
plot( PDF2(:,2), PDF2(:,1), 'r' );

title( 'Histograms and Pdf''s ' );
xlabel( 'Random Variable (Unitless)' );
ylabel( 'Probability Density' );
grid
legend( 'Histogram', 'Mean PDF');
print -dpng Figure2.png

%% 
% for row 3
% Plot the histogram estimate of the Histogram
figure(3);
bar( HistogramPDF3(:,2), HistogramPDF3(:,1) );
hold on

% Plot the pdf's using the Gaussian distribution of 90% CI
plot( PDF3(:,2), PDF3(:,1), 'r' );

title( 'Histograms and Pdf''s ' );
xlabel( 'Random Variable (Unitless)' );
ylabel( 'Probability Density' );
grid
legend( 'Histogram', 'Mean PDF');
print -dpng Figure3.png

%% 
% for row 4
% Plot the histogram estimate of the Histogram
figure(4);
bar( HistogramPDF4(:,2), HistogramPDF4(:,1) );
hold on

% Plot the pdf's using the Gaussian distribution of 90% CI
plot( PDF4(:,2), PDF4(:,1), 'r' );

title( 'Histograms and Pdf''s ' );
xlabel( 'Random Variable (Unitless)' );
ylabel( 'Probability Density' );
grid
legend( 'Histogram', 'Mean PDF');
print -dpng Figure4.png

%% 
% for row 5
% Plot the histogram estimate of the Histogram
figure(5);
bar( HistogramPDF5(:,2), HistogramPDF5(:,1) );
hold on

% Plot the pdf's using the Gaussian distribution of 90% CI
plot( PDF5(:,2), PDF5(:,1), 'r' );

title( 'Histograms and Pdf''s ' );
xlabel( 'Random Variable (Unitless)' );
ylabel( 'Probability Density' );
grid
legend( 'Histogram', 'Mean PDF');
print -dpng Figure5.png
