# RMTWasserstein
Code for the paper Random  Matrix-Improved  Estimation  of  the  Wasserstein  Distance between  two  Centered  Gaussian  Distribution.

The code contains the following files:

The RWTWassDist which is a function which takes as arguments the two samples X and Y from the 2 classes and return the improved estimate est and the classical estimate esthat.

The script test_wasserstein which allows to reproduce the table 1 provided in the paper. it compares the classical estimate to the proposed one.

The script CompareEst is the application to the covariance matrix estimation through the gradient descent algorithm.

The others are utilities for the gradient descent algorithms.

Feel free to contact the authors at tiomoko_malik@yahoo.fr for further details.
