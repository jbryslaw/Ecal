# Ecal
Calibration of Electromagnetic Calorimeters Using Machine Learning Methods

+ generate_data.C
Generates test and training data assuming signal is of the form:
sig(x) = Landau_0(mu,sigma,x) + Landau_1(Landau_0,x)

and background is of the form:
Bb(x) = Gauss(x)+exp(x)

with a random square cutout at low x to simulate zero-suppression.

Default number of trials to be divided among training and testing: 200000