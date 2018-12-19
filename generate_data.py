# save histogram bin content and fit parameters
#   to txt or binary file instead of saving histogram objects

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import warnings; warnings.filterwarnings(action='once')

#needed for defining landau distributions
import pylandau as pylan #pip install pylandau


#define a linspace for the ADC distributions
x_space = np.linspace(0,150,150)
mpv, sig, amp = 30,5,100
#test_landau = pylan.landau(x_space,mpv = mpv, eta=sig, A=amp)
#test_landau = [pylan.get_landau(d_ijk,mpv=mpv,eta=sig,A=amp) for d_ijk in x_space]
#test_landau = [pylan.get_landau(d_ijk,mpv=30,eta=5,A=10) for d_ijk in range(0,150)]
#test_landau = [pylan.landau(np.arange(d_ijk,d_ijk+1.,1),mpv=30,eta=5,A=10) for d_ijk in range(0,150)]
test_landau = [pylan.landau(np.array([d_ijk+0.0]),mpv=30,eta=5,A=10) for d_ijk in range(0,150)]
#plt.plot(x_space,test_landau)
#plt.hist(test_landau)
print([d_ijk for d_ijk in x_space])

# Energy Deposition:
# 9 parameter function Gaussian + Landau0 + Landa1(Landa0) + exp
def eval_E_deposit_func(x,par):
    x = x*1.
    if len(par) != 9:
        return 0.

    for ijk in par:
        #first Landau:
        #la0 = par[3]*pylan.get_landau(x,par[1],par[2],A=1) # get_landau function doesn't seem to return the proper mpv
        #use the "full" distribution with one x entry instead
        la0 = par[3]*pylan.landau(np.array([x]),mpv = par[1], eta = par[2])
        la1 = (1.-par[3])*pylan.landau(np.array([x]),mpv = 2*par[1]+1.4*par[2],eta = 2.*par[2])

        bg  = par[4]*np.exp(-1.*par[5]*x)+par[6]*np.exp(-0.5*((x-par[7])/par[8])**2)

        total = par[0]*(la0[0]+la1[0])+bg
    return total

par = [0,1,2,3,4,5,6,7,8]

par = [1.6e4,
       18.92,
       2.5,
       1-0.09,
       351,
       0.1,
       1e4,
       -20,
       9]

x_val = [ i+0.5 for i in range(0,150)]
y_val = [ eval_E_deposit_func(x,par) for x in x_val]
print(eval_E_deposit_func(23,par))
# print(x_val)
# print(y_val)
plt.plot(x_val,y_val)
