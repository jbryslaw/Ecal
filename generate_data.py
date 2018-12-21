# save histogram bin content and fit parameters
#   to txt or binary file instead of saving histogram objects

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import warnings; warnings.filterwarnings(action='once')
from scipy import stats as scipystats

#################
# Switches
b_plot = False #True
b_savetxt = False
b_saveplk = True
#################

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
        return np.floor(x/self.binwidth).astype(int)

#class th1d_hist:
#####################################################################

#define a linspace for the ADC distributions
x_val = [ i+0.5 for i in range(0,150)]

#####################################################################
# Energy Deposition:
# 9 parameter function Gaussian + Landau0 + Landa1(Landa0) + exp
def eval_E_deposit_func(x,par):
    x = x*1.
    if len(par) != 9:
        return 0.
    
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

np.random.seed(0)

#Generate 100000 training and test histograms
x_max = 150.
i_Nsamples = 100000
i_N_itr = 200000;

#limit for testing
i_Nsamples = 100000
i_Nsamples = 10
i_N_itr = 10

histos = []


# save fit parameters for fit
# and histogram bins for data
l_columns = [f"p{i}"  for i in range(0,9)]
l_bins_columns = [f"b{i}"  for i in range(0,150)]

l_columns.extend(l_bins_columns)

#df_total = pd.DataFrame(columns=l_columns)
for ijk in range(0,i_N_itr):
    #get a random set of parameters
    par=[np.random.rand()*(par_hi[jkl]-par_low[jkl])+par_low[jkl] for jkl in range(0,len(par_low))]
    a_par = np.array(par)
            
    # randomly sample function
    h1 = th1d_hist(150,0.,150.)
    for jkl in range(0,i_Nsamples):
        xsample = np.random.rand()
        xsample *= x_max
        weight  = eval_E_deposit_func(xsample,par)
        h1.Fill(xsample,weight)
    #for jkl in range(0,i_Nsamples):

    #fill random for testing
    # for jkl in range(0,10000):
    #     h1.Fill(np.floor(np.random.rand()*150).astype(int),np.random.rand() )

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
        
    histos.append(h1)

    #save fit parameters and bin content to dataframe
    parbin_out = par
    ls_content = [h1.GetBinContent(ijk) for ijk in range(0,150)]
    parbin_out.extend(ls_content)
    if ijk == 0:    
        #df_total = pd.DataFrame(a_par,columns=l_columns)
        df_total = pd.DataFrame([parbin_out],columns = l_columns)
    else:
        df0 = pd.DataFrame([parbin_out],columns=l_columns)
        df_total =  pd.concat([df_total,df0], ignore_index=True)


#for ijk in range(0,i_N_itr)
print(df_total)

# plot all histograms
if b_plot:
    for ijk in histos:
        plt.figure()
        plt.semilogy(x_val,ijk.get_hist_array())
#if b_plot:

#save as pickle
if b_saveplk:
    df_total.to_pickle("all_data.plk")

