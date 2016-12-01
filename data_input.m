% Given data in the assignment
% Just run the matlab-file to load the data to another script

% Covariance
H = [4.01 -1.19 0.60 0.74 -0.21;
     -1.19 1.12 0.21 0.54 0.55;
      0.60 0.21 3.04 0.77 0.29;
      0.74 0.54  0.77 3.74 -1.04 
     -0.21 0.55  0.29 -1.04  3.80]*1e-2;

% Average return
r = [13.0, 5.3, 10.5, 5.0, 12.6]*1e-2;

% Standard deviation of each asset
stdev = sqrt(diag(H));