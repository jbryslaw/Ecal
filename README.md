# Ecal
Calibration of Electromagnetic Calorimeters Using Machine Learning Methods

0. generate_data.C:  
Generates test and training data with Monte Carlo assuming signal is of the form:
sig(x) = Landau_0(mu,sigma,x) + Landau_1(Landau_0,x)

and background is of the form:\
Bb(x) = Gauss(x)+exp(x)

with a random square cutout at low x to simulate zero-suppression.\

Default number of trials to be divided among training and testing: 200000\
Output is in the form of a .root file containing generated histograms along with original generating functions.

1. get_features.C:\
Finds features of inputed histograms. Secifically, find local maxima, minima and locations where dydx=0 and d2ydx2=0.c

2. find_mpv_KNN.R\
Uses KNN to find MIP MPV from histogram features.