#include <iostream>
#include "Riostream.h"
//#include <unistd.h>
#include <cstdlib>
//#include <pthread.h>
// #include <cstdlib>
// #include <pthread.h>
// //#include <iomanip>
#include <string>
#include <vector>
// //#include <algorithm>
// //#include <cstdlib>
// //#include "gsl/gsl_rng.h"
#include <math.h>
#include <set>
#include <map>

#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStreamerElement.h>
#include <TStyle.h>
#include "TSystemDirectory.h"

#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGaxis.h"
#include "TSpectrum.h"
// #include "TLorentzVector.h"
// #include "TVector3.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TChain.h"
#include "TLegend.h"
#include "TFractionFitter.h"
#include "TVirtualFitter.h"
#include "TCut.h"
#include "TObject.h"
#include "TLine.h"
#include "TSpline.h"
#include "TPaveText.h"
#include "TFormula.h"

#include "Fit/FitResult.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "TRandom3.h"
#include "TDatime.h"

#include "Math/MinimizerOptions.h"

using namespace std;

void generate_data()
{
  // Output Histograms to Screen?
  bool b_draw = true;

  // Create Files
  TFile * tf_train_out;
  TFile * tf_test_out;

  if(!b_draw)
    {
      tf_train_out = new TFile("training_data.root","RECREATE"); 
      tf_test_out  = new TFile("testing_data.root","RECREATE");  
    }//if(!b_draw)

  //Signal: f(x) = a*Landau(mu_0,sigma_0) + (1-a)Landau(2*mu_0+1.4*sigma_0,2*sigma_0)
  const char *la0  = "[3]*TMath::Landau(x,[1],[2],0)";
  const char *la1  = "(1-[3])*TMath::Landau(x,2*[1]+1.4*[2],2.0*[2],0)";

  //Background: f(x) = a*Gaus(mu_1,sigma_1) + b*exp(tau)
  const char *c_bg = "[4]*exp(-[5]*x)+[6]*exp(-0.5*((x-[7])/[8])**2)";

  // zero suppression: random -step function at low x

  // Set Parameter Ranges from Real Data
  double a_d_par_low[9] =
    {
      1.6e4,
      18.92,
      2.5,
      1-0.09,
      351,
      0.1,
      1e4,
      -20,
      9
    };//double a_d_par_low[9] =

  double a_d_par_hi[9] =
    {
      8.8e5,
      24.13,
      3.6,
      1-0.05,
      1.742e4,
      0.23,
      2e5,
      10,
      15      
    };//double a_d_par_hi[9] =

  double a_d_par[9];

  TFormula * tform_s_bg = new TFormula("tform_s_bg", Form("[0]*(%s+%s)+%s",la0,la1,c_bg));

  TRandom3 * tr3 = new TRandom3();
  tr3 -> SetSeed(1);

  // Count the number of saved histograms
  TH1I * h_histos_train = new TH1I("h_histos_train","Number of Histograms",1,0,1);
  TH1I * h_histos_test  = new TH1I("h_histos_test","Number of Histograms",1,0,1);

  // Generate 100000 training and test histograms
  int i_N_itr = 200000;
  for(int ijk = 0; ijk < i_N_itr; ijk ++)
    {
      // Get a Function from Set Formula
      TF1 * f1 = new TF1(Form("f1_%d",ijk),"tform_s_bg",0,150);

      // Set random parameters within ranges
      for(int jkl = 0; jkl < 9; jkl ++)
	{
	  a_d_par[jkl] = (tr3 -> Rndm())*(a_d_par_hi[jkl]-a_d_par_low[jkl])+a_d_par_low[jkl];
	  f1 -> SetParameter(jkl,a_d_par[jkl]);
	}//for(int jkl = 0; jkl < 9; jkl ++)  

      // Save functon integral for later
      double d_range0 = 0.0;
      double d_range1 = 150.0;
      double d_func_int = f1->Integral(d_range0,d_range1);
      if(d_func_int == 0.0) continue;

      // Fill Histogram with random samples from generator function
      TH1D * h_energy = new TH1D(Form("h_energy_%d",ijk),"Energy Deposition",150,0,150);
      h_energy -> Sumw2();
      h_energy -> FillRandom(Form("f1_%d",ijk),100000);

      // Normalize histogram to generator function
      int iBin0 = h_energy -> FindFixBin(d_range0);
      int iBin1 = h_energy -> FindFixBin(d_range1);
      double d_hist_int = h_energy -> Integral(iBin0,iBin1);
      h_energy -> Scale(d_func_int/d_hist_int);

      // zero suppression
      // Random Rectangular Cut Out from 0 to 20
      double d_zs_range    = (tr3->Rndm()*20.0);
      double d_zs_strength = 10*tr3->Rndm();

      double d_bin0_content = h_energy -> GetBinContent(1);
      int iBin_zs_range = h_energy -> FindFixBin(d_zs_range);
      for(int jkl = 1; jkl < iBin_zs_range; jkl++) h_energy -> SetBinContent(jkl,d_bin0_content*pow(10,-d_zs_strength));

      //count training and test histograms
      if(!b_draw)
	{
	  if((2*ijk) < i_N_itr)
	    {
	      tf_train_out -> cd();
	      h_histos_train -> SetBinContent(1,h_histos_train->GetBinContent(1)+1);		
	    }//if((2*ijk) < i_N_itr)
	  else
	    {
	      tf_test_out -> cd();
	      h_histos_test -> SetBinContent(1,h_histos_test->GetBinContent(1)+1);	  
	    }//else ((2*ijk) < i_N_itr)
	  //write to file
	  f1 -> Write();
	  h_energy -> Write();
	}//if(!b_draw)
      
      // Draw
      if(b_draw)
	{
	  new TCanvas();
	  gPad -> SetLogy();
	  h_energy->Draw("HIST");
	  f1 -> Draw("same");
	  if(ijk > 2) return; // don't draw too many histograms to screen
	}//if(b_draw)
    }//for(int ijk = 0; ijk < i_N_itr; ijk ++)

  tf_train_out -> cd();
  h_histos_train -> Write("h_histos");
  tf_test_out -> cd();
  h_histos_test -> Write("h_histos");

  return;
}//void generate_data()


//is only read by compiler
#ifndef __CINT__

int main(int argc, char * argv[]){
  gSystem->Load("libTree");
  TApplication App("Analysis", &argc, argv);
  
  generate_data();
  
  //App.Terminate();
  cout << "Main: program exiting." << endl;
  pthread_exit(NULL);

  return 0;
}
#endif
