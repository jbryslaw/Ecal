import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import warnings; warnings.filterwarnings(action='once')
from scipy import stats as scipystats
import pylab

#########################
#  switches
b_draw             = True #False
b_draw_derivatives = True #False #True #False
#########################

#####################################################################
# define histogram like object:
class th1d_hist:
    def __init__(self,i_nBins,x0,x1,title="title"):
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

    def GetBinCenter(self,i_bin):
        return self.x_val[i_bin]

    def GetNbins(self):
        return len(self.a_hist)

    def GetBinWidth(self):
        return self.binwidth

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
            self.SetBinContent(ix,w)
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
i_N_rows = df_total.shape[0]

i_max_rows = 2
i_start = 0

for ijk in range(i_start, i_N_rows):
    print(" ------ ",f'hist_{ijk}')
    if(ijk >= i_max_rows):
        break

    # get histogram bins
    this_hist = df_total.iloc[ijk,9:159]
    h_energy = th1d_hist(150,0.,150.,f'hist_{ijk}')
    h_energy.read(this_hist)
    h_energy.set_poisson_err(100000)

    # Draw
    if b_draw:
        plt.figure(figsize=(16,10),dpi=80)
        ax = plt.subplot(111)
        ax.set_yscale("log")
        h_energy.Draw()
    #if b_draw:

    #--------------------------------------------
    # Local Minima and Maxima
    v_iBin_minima = list()
    v_iBin_maxima = list()
    
    # loop over bins
    for iBin in range(0,h_energy.GetNbins()):
        d_center  = h_energy.GetBinCenter(iBin)
        d_content = h_energy.GetBinContent(iBin)

        #Define Max, Min as 5 bins lower, higher on either side
        #  so statistical fluctuations aren't set as max,min
        a5_pre_content  = [0.,0.,0.,0.,0.]
        a5_post_content = [0.,0.,0.,0.,0.]

        b_min = True
        b_max = True
        i_num_check = 5
        
        for jkl in range(0,i_num_check):
            #just check that bin isn't near the very edges of the histogram
            if (iBin-1-jkl) > 0:
                a5_pre_content[jkl] = h_energy.GetBinContent(iBin-1-jkl)
            else:
                b_min = False
                b_max = False
            if (iBin+1+jkl) < h_energy.GetNbins():
                a5_post_content[jkl] = h_energy.GetBinContent(iBin+1+jkl)
            else:
                b_min = False
                b_max = False

            if(a5_pre_content[jkl] < d_content): b_min = False
            if(a5_pre_content[jkl] > d_content): b_max = False
            if(a5_post_content[jkl] < d_content): b_min = False
            if(a5_post_content[jkl] > d_content): b_max = False

            #if the bins right nextdoor are the same, let it go
            if(jkl==0): continue
            if(a5_pre_content[jkl]  == d_content): b_max = False
            if(a5_post_content[jkl] == d_content): b_min = False
        #for jkl in range(0,i_num_check):

        #add maxima, minima to vectors
        if(b_max):
            v_iBin_maxima.append(iBin)
                
        if(b_min):
            v_iBin_minima.append(iBin)
    #for iBin in range(0,h_energy.GetNbins()):

    #draw arrows at minima nad maxima
    print("Minima:")
    for jkl in range(0,len(v_iBin_minima)):
        print(v_iBin_minima[jkl])
        d_center  = h_energy.GetBinCenter(v_iBin_minima[jkl])
        d_content = h_energy.GetBinContent(v_iBin_minima[jkl])
        #pylab.arrow(d_center,d_content-(d_content*0.1),0.0,(d_content*0.1),fc="b",ec="b",head_width=0.1, head_length=d_content )
        plt.annotate('',xy=(d_center,d_content),xytext=(d_center,d_content-(d_content*0.1)),arrowprops=dict(facecolor='blue',shrink=0.))

    print("Maxima:")
    for jkl in range(0,len(v_iBin_maxima)):
        print(v_iBin_maxima[jkl])
        d_center  = h_energy.GetBinCenter(v_iBin_maxima[jkl])
        d_content = h_energy.GetBinContent(v_iBin_maxima[jkl])
        #pylab.arrow(d_center,d_content-(d_content*0.1),0.0,(d_content*0.1),fc="b",ec="b",head_width=0.1, head_length=d_content )
        plt.annotate('',xy=(d_center,d_content),xytext=(d_center,d_content-(d_content*0.1)),arrowprops=dict(facecolor='red',shrink=0.))

    # END Local Minima and Maxima
    #--------------------------------------------

    #--------------------------------------------
    #       Local Derivatives

    #first find all derivatives and save to a histogram
    h_dydx = th1d_hist(150,0,150)
    # loop over bins
    for iBin in range(0,h_energy.GetNbins()-1):
        # Skip edges of histogram
        if((iBin < 1) or (iBin > h_energy.GetNbins())):
            continue

        dx        = h_energy.GetBinWidth()
        if(dx==0):
            continue
        d_center  = h_energy.GetBinCenter(iBin)
        d_content = h_energy.GetBinContent(iBin)

        d_pre_content  = h_energy.GetBinContent(iBin-1)
        d_next_content = h_energy.GetBinContent(iBin+1)

        dydx = (d_next_content - d_pre_content) /dx
        h_dydx.SetBinContent(iBin,dydx)
    #for iBin in range(0,h_energy.GetNbins()):

    if(b_draw_derivatives):
        plt.figure(figsize=(16,10),dpi=80)
        h_dydx.Draw()

    #find dy/dx = 0
    v_iBin_dxdy0 = list()
    for iBin in range(1,h_dydx.GetNbins()-1):
        #start assuming dydx != 0
        b_dydx0 = False

        d_center  = h_dydx.GetBinCenter(iBin)
        d_content = h_dydx.GetBinContent(iBin)
        d_pre_content  = h_dydx.GetBinContent(iBin-1)
        d_next_content = h_dydx.GetBinContent(iBin+1)

        if((d_content < 0) and (d_pre_content > 0)): b_dydx0 = True
        if((d_content > 0) and (d_pre_content < 0)): b_dydx0 = True

        if(not b_dydx0): continue
        v_iBin_dxdy0.append(iBin)

        #Draw arrows
        # assuming h_energy and h_dydx have
        #     the same number of bins
        plt.annotate('',xy=(d_center,d_content),xytext=(d_center,d_content-(d_content*0.1)),arrowprops=dict(facecolor='green',shrink=0.))
        
    #   END Local Derivatives
    #--------------------------------------------

    #--------------------------------------------
    #       Local 2nd Derivatives
    #first find all 2nd derivatives and save to a histogram
    h_d2ydx2 = th1d_hist(150,0,150)

    #loop over bins in the 1st derivative histogram
    for iBin in range(0,h_dydx.GetNbins()-1):
        # Skip edges of histogram
        if((iBin < 1) or (iBin > h_dydx.GetNbins())):
            continue

        dx        = h_dydx.GetBinWidth()
        if(dx==0):
            continue
        d_center  = h_dydx.GetBinCenter(iBin)
        d_content = h_dydx.GetBinContent(iBin)

        d_pre_content  = h_dydx.GetBinContent(iBin-1)
        d_next_content = h_dydx.GetBinContent(iBin+1)

        dydx = (d_next_content - d_pre_content) /dx
        h_d2ydx2.SetBinContent(iBin,dydx)
    #for iBin in range(0,h_dydx.GetNbins()):

    if(b_draw_derivatives):
        plt.figure(figsize=(16,10),dpi=80)
        h_d2ydx2.Draw()

    #find d2y/dx2 = 0
    v_iBin_d2ydx2_0 = list()
    for iBin in range(1,h_d2ydx2.GetNbins()-1):
        # start assuming d2y/dx2 != 0
        b_d2ydx2_0 = False
        d_center  = h_d2ydx2.GetBinCenter(iBin)
        d_content = h_d2ydx2.GetBinContent(iBin)

        d_pre_content   = h_d2ydx2.GetBinContent(iBin-1)
        d_next_content  = h_d2ydx2.GetBinContent(iBin-1)

        #Appears to work best if this bins is compared to the next bin
        if((d_content < 0) and (d_pre_content > 0)): b_d2ydx2_0 = True
        if((d_content > 0) and (d_pre_content < 0)): b_d2ydx2_0 = True

        if(not b_d2ydx2_0): continue
        v_iBin_d2ydx2_0.append(iBin)

        #Draw arrows
        # assuming h_energy and h_dydx have
        #     the same number of bins
        plt.annotate('',xy=(d_center,d_content),xytext=(d_center,d_content-(d_content*0.1)),arrowprops=dict(facecolor='black',shrink=0.))

        
    #   END Local 2nd Derivatives
    #--------------------------------------------



#for ijk in range(0, i_N_rows):
