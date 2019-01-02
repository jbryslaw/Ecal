import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import warnings; warnings.filterwarnings(action='once')
from scipy import stats as scipystats

#########################
#  switches
b_draw = True
#########################

#####################################################################
# define histogram like object:
class th1d_hist:
    def __init__(self,i_nBins,x0,x1,title):
        self.i_nBins = i_nBins
        self.x0 = x0
        self.x1 = x1
        self.title = title

        self.disp = x1-x0
        self.binwidth = -9999.0
        if i_nBins != 0:
            self.binwidth = self.disp/(i_nBins*1.)

        #self.a_hist = np.array([x0+(i*self.binwidth) for i in range(0,i_nBins)])
        self.a_hist = np.zeros(i_nBins)
        self.a_err  = np.zeros(i_nBins)
        self.x_val  = np.arange(self.x0,self.x1,self.binwidth)

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

    #set histogram poisson errors
    def set_poisson_err(self,d_norm):
        a_ones = np.ones(self.i_nBins)
        self.a_err = np.zeros(self.i_nBins)
        a_dsq = np.zeros(self.i_nBins)
        a_dsq = np.sqrt(self.a_hist*d_norm)
        self.a_err = np.true_divide(a_ones,a_dsq,out=np.zeros_like(a_ones),where=(a_dsq!=0.))

    # read in histogram values directly
    # assumes source has same number of entries as histogram has bins
    def read(self,v_in):
        ix = 0
        for w in v_in:            
            h1.SetBinContent(ix,w)
            ix+=1

    def get_xarray(self):
        return self.x_val

    # Draw
    def Draw(self):
        half_binwidth = self.binwidth/2.
        x_draw = np.arange(self.x0+half_binwidth,self.x1+half_binwidth,self.binwidth)
        plt.errorbar(x_draw,self.a_hist,yerr=self.a_err,fmt='o',  markersize=2, color='k', label = self.title)
        
#class th1d_hist:
#####################################################################


#define a linspace for the ADC distributions
x_val = [ i+0.5 for i in range(0,150)]
# read in the data file
df_total = pd.read_pickle("all_data.plk")

print(" shape all data: ")
print(df_total.shape)

# Number of historgrams to read
# int i_N_histos = 1000

# for ijk in range(0,i_N_histos):

# get histogram bins
this_hist = df_total.iloc[0,9:159]
h1 = th1d_hist(150,0.,150.,"hist0")
h1.read(this_hist)
h1.set_poisson_err(10000)

# Draw

h1.Draw()

# print(this_hist.shape)
# print(this_hist)

# h1.Fill(x_val,this_hist)
    

# plt.figure(figsize=(16,10),dpi=80)
# #plt.scatter(this_hist.index,this_hist)
# #plt.errorbar(x_val,this_hist,
