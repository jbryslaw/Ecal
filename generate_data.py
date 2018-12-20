# save histogram bin content and fit parameters
#   to txt or binary file instead of saving histogram objects

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import warnings; warnings.filterwarnings(action='once')
from scipy import stats as scipystats

#needed for defining landau distributions --> use scipy.stats.moyal instead
#import pylandau as pylan #pip install pylandau

#####################################################################
# define histogram like object:
class th1d_hist:
    def __init__(self,i_nBins,x0,x1):
        self.i_nBins = i_nBins
        self.x0 = x0
        self.x1 = x1

        self.disp = x1-x0
        self.binwidth = -9999.0
        if i_nBins != 0:
            self.binwidth = self.disp/(i_nBins*1.)

        #self.a_hist = np.array([x0+(i*self.binwidth) for i in range(0,i_nBins)])
        self.a_hist = np.zeros(i_nBins)

    def get_hist_array(self):
        return self.a_hist

    def Fill(self,x=0.,w=1.):
        i_bin = np.floor(x/self.binwidth).astype(int)
        self.a_hist[i_bin] += w

    def Scale(self,scale):
        self.a_hist *= scale

    def SetBinContent(self,i_bin,d_content):
        self.a_hist[i_bin] = d_content

    def GetBinContent(self,i_bin):
        return self.a_hist[i_bin]

    def FindBin(self,x):
        print(" - ", x," - ",self.binwidth," - ", np.floor(x/self.binwidth), " - ", np.floor(x/self.binwidth).astype(int) )
        return np.floor(x/self.binwidth).astype(int)

#class th1d_hist:
#####################################################################

#define a linspace for the ADC distributions
x_space = np.linspace(0,150,150)
mpv, sig, amp = 30,5,100
#test_landau = pylan.landau(x_space,mpv = mpv, eta=sig, A=amp)
#test_landau = [pylan.get_landau(d_ijk,mpv=mpv,eta=sig,A=amp) for d_ijk in x_space]
#test_landau = [pylan.get_landau(d_ijk,mpv=30,eta=5,A=10) for d_ijk in range(0,150)]
#test_landau = [pylan.landau(np.arange(d_ijk,d_ijk+1.,1),mpv=30,eta=5,A=10) for d_ijk in range(0,150)]
#test_landau = [pylan.landau(np.array([d_ijk+0.0]),mpv=30,eta=5,A=10) for d_ijk in range(0,150)]
#plt.plot(x_space,test_landau)
#plt.hist(test_landau)
#print([d_ijk for d_ijk in x_space])

#####################################################################
# Energy Deposition:
# 9 parameter function Gaussian + Landau0 + Landa1(Landa0) + exp
def eval_E_deposit_func(x,par):
    x = x*1.
    if len(par) != 9:
        return 0.

    #first Landau:
    #la0 = par[3]*pylan.get_landau(x,par[1],par[2],A=1) # get_landau function doesn't seem to return the proper mpv
    #use the "full" distribution with one x entry instead
    # la0 = par[3]*pylan.landau(np.array([x]),mpv = par[1], eta = par[2])
    # la1 = (1.-par[3])*pylan.landau(np.array([x]),mpv = 2*par[1]+1.4*par[2],eta = 2.*par[2])
    
    # can use moyal distributrion from scipy.stats to approximate landau, it is much faster
    la0 = par[3]*scipystats.moyal.pdf(x,loc=par[1],scale=par[2])
    la1 = (1.-par[3])*scipystats.moyal.pdf(x,loc=(2.*par[1]+1.4*par[2]), scale= (2.*par[2]))
    bg  = par[4]*np.exp(-1.*par[5]*x)+par[6]*np.exp(-0.5*((x-par[7])/par[8])**2)

    total = par[0]*(la0+la1)+bg
    return total
#def eval_E_deposit_func(x,par):
#####################################################################
par_low = [1.6e4,
           18.92,
           2.5,
           1-0.09,
           351,
           0.1,
           1e4,
           -20,
           9]

par_hi = [8.8e5,
          24.13,
          3.6,
          1-0.05,
          1.742e4,
          0.23,
          2e5,
          10,
          15]


x_val = [ i+0.5 for i in range(0,150)]
#y_val = [ eval_E_deposit_func(x,par) for x in x_val]
#print(eval_E_deposit_func(23,par))
# print(x_val)
# print(y_val)
b_plot = False
if b_plot:
    plt.semilogy(x_val,y_val)


np.random.seed(0)
    
# plt.semilogy(x_val,h1.get_hist_array())
#plt.plot(x_val,h1.get_hist_array())

#Generate 100000 training and test histograms
x_max = 150.
i_Nsamples = 100000
i_N_itr = 200000;

#limit for testing
i_Nsamples = 100000
i_N_itr = 10

histos = []

# save fit parameters for fit
# and histogram bins for data
for ijk in range(0,i_N_itr):
    #get a random set of parameters
    par=[np.random.rand()*(par_hi[jkl]-par_low[jkl])+par_low[jkl] for jkl in range(0,len(par_low))]

    # randomly sample function
    h1 = th1d_hist(150,0.,150.)
    for jkl in range(0,i_Nsamples):
        xsample = np.random.rand()
        xsample *= x_max
        weight  = eval_E_deposit_func(xsample,par)
        h1.Fill(xsample,weight)
    #for jkl in range(0,i_Nsamples):

    if i_Nsamples != 0. :
        h1.Scale(1./(i_Nsamples*1.))

    # zero suppression
    # Random Rectangular Cut Out from 0 to 20
    zs_range    = np.random.rand()*20
    zs_strength = 10*np.random.rand()
    d_old_bincontent = h1.GetBinContent(1)
    zs_rightbin = h1.FindBin(zs_range)
    #print("zs_rightbin: ",zs_rightbin, " zs_strength ",zs_strength)

    for jkl in range(0,zs_rightbin):
        h1.SetBinContent(jkl,d_old_bincontent*np.power(10,-1*zs_strength))
        print(" - ",jkl," ",d_old_bincontent*np.power(10,-1*zs_strength))

    histos.append(h1)
#for ijk in range(0,i_N_itr)

for ijk in histos:
    plt.figure()
    plt.semilogy(x_val,ijk.get_hist_array())
