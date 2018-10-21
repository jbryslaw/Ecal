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

//--- 
//need this stuff to compile in linux:
#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStreamerElement.h>
#include <TStyle.h>
#include "TSystemDirectory.h"
//---

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
  TFile * tf_train_out = new TFile("training_data.root","RECREATE");
  TFile * tf_test_out  = new TFile("testing_data.root","RECREATE");
   
  //Signal: f(x) = a*Landau(mu_0,sigma_0) + (1-a)Landau(2*mu_0+1.4*sigma_0,2*sigma_0)
  //        15<mu_0<25 3<sigma_0<6
  //Background: f(x) = a*Gaus(mu_1,sigma_1) + b*exp(tau)
  // zero suppression: random -step function at low x

  const char *la0  = "[3]*TMath::Landau(x,[1],[2],0)";
  const char *la1  = "(1-[3])*TMath::Landau(x,2*[1]+1.4*[2],2.0*[2],0)";
  const char *c_bg = "[4]*exp(-[5]*x)+[6]*exp(-0.5*((x-[7])/[8])**2)";

  //ranges
  double d_par0_0 = 1.6e4;
  double d_par0_1 = 8.8e5;
  double d_par1_0 = 18.92;
  double d_par1_1 = 24.13;
  double d_par2_0 = 2.5;
  double d_par2_1 = 3.6;
  double d_par3_0 = 1-0.09;
  double d_par3_1 = 1-0.05;
  double d_par4_0 = 351;
  double d_par4_1 = 1.742e4;
  double d_par5_0 = 0.1;
  double d_par5_1 = 0.23;
  double d_par6_0 = 1e4;
  double d_par6_1 = 2e5;
  double d_par7_0 = -20;
  double d_par7_1 = 10;
  double d_par8_0 = 9;
  double d_par8_1 = 15;

  TFormula * tform_s_bg = new TFormula("tform_s_bg", Form("[0]*(%s+%s)+%s",la0,la1,c_bg));

  TRandom3 * tr3 = new TRandom3();

  //tr3 -> Rndm();
  vector<TH1D*> v_energy_histos;

  TH1I * h_histos_train = new TH1I("h_histos_train","Number of Histograms",1,0,1);
  TH1I * h_histos_test  = new TH1I("h_histos_test","Number of Histograms",1,0,1);
  
  int i_N_itr = 200000;
  for(int ijk = 0; ijk < i_N_itr; ijk ++)
    {
      //cout<<" "<<d_par0_0<<" "<<d_par0_1<<" "<<(tr3 -> Rndm())*(d_par0_1-d_par0_0)+d_par0_0<<endl;
      double d_par0 = (tr3 -> Rndm())*(d_par0_1-d_par0_0)+d_par0_0;
      double d_par1 = (tr3 -> Rndm())*(d_par1_1-d_par1_0)+d_par1_0;
      double d_par2 = (tr3 -> Rndm())*(d_par2_1-d_par2_0)+d_par2_0;
      double d_par3 = (tr3 -> Rndm())*(d_par3_1-d_par3_0)+d_par3_0;
      double d_par4 = (tr3 -> Rndm())*(d_par4_1-d_par4_0)+d_par4_0;
      double d_par5 = (tr3 -> Rndm())*(d_par5_1-d_par5_0)+d_par5_0;
      double d_par6 = (tr3 -> Rndm())*(d_par6_1-d_par6_0)+d_par6_0;
      double d_par7 = (tr3 -> Rndm())*(d_par7_1-d_par7_0)+d_par7_0;
      double d_par8 = (tr3 -> Rndm())*(d_par8_1-d_par8_0)+d_par8_0;

      
      TF1 * f1 = new TF1(Form("f1_%d",ijk),"tform_s_bg",0,150);
      f1 -> SetParameter(0,d_par0);
      f1 -> SetParameter(1,d_par1);//25.0);
      f1 -> SetParameter(2,d_par2);//2.60);
      f1 -> SetParameter(3,d_par3);//0.8);
      f1 -> SetParameter(4,d_par4);//1.0);
      f1 -> SetParameter(5,d_par5);//1.0);
      f1 -> SetParameter(6,d_par6);//0.06);
      f1 -> SetParameter(7,d_par7);//8000.0);
      f1 -> SetParameter(8,d_par8);//0.02);

      double d_range0 = 0.0;
      double d_range1 = 150.0;
      double d_func_int = f1->Integral(d_range0,d_range1);
      if(d_func_int == 0.0) continue;

      TH1D * h_energy = new TH1D(Form("h_energy_%d",ijk),"Energy Deposition",150,0,150);
      h_energy -> Sumw2();
      h_energy -> FillRandom(Form("f1_%d",ijk),100000);

      int iBin0 = h_energy -> FindFixBin(d_range0);
      int iBin1 = h_energy -> FindFixBin(d_range1);
      double d_hist_int = h_energy -> Integral(iBin0,iBin1);
      h_energy -> Scale(d_func_int/d_hist_int);

      // zero suppression

      double d_zs_range    = (tr3->Rndm()*20.0);
      double d_zs_strength = 10*tr3->Rndm();

      double d_bin0_content = h_energy -> GetBinContent(1);
      int iBin_zs_range = h_energy -> FindFixBin(d_zs_range);
      for(int jkl = 1; jkl < iBin_zs_range; jkl++) h_energy -> SetBinContent(jkl,d_bin0_content*pow(10,-d_zs_strength));

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

      f1 -> Write();
      h_energy -> Write();

      // Draw
      //      new TCanvas(); gPad -> SetLogy();  h_energy->Draw("HIST"); f1 -> Draw("same");
    }//for(int ijk = 0; ijk < i_N_itr; ijk ++)

  tf_train_out -> cd();
  h_histos_train -> Write("h_histos");
  tf_test_out -> cd();
  h_histos_test -> Write("h_histos");

  // TF1 * f1 = new TF1("f1","tform_s_bg",0,150);
  // f1 -> SetParameter(0,4.0e5);
  // f1 -> SetParameter(1,25.0);
  // f1 -> SetParameter(2,2.60);
  // f1 -> SetParameter(3,0.8);
  // f1 -> SetParameter(4,1.0);
  // f1 -> SetParameter(5,1.0);
  // f1 -> SetParameter(6,0.06);
  // f1 -> SetParameter(7,8000.0);
  // f1 -> SetParameter(8,0.02);

  //  f1 -> SetParameter(0,0.0);
  //f1 -> SetParameter(6,0.0);

  //with exp
  //0: 1.6e4 : 1.7e4 1.78e4 : 3.39e5 :5.4e5 : 7.4e5 
  //1: 13.31 : 18.92 : 19.51 : 21.1 : 24.13
  //2: 2.5 : 2.781 : 2.939 : 3.00 : 3.578
  //3: 1-[] 0.0562 :0.06692 :  0.072 : 0.076 : 0.09	   
  //4: 351 : 439 :732 : 6525 : 1.742e4		  
  //5: 0.011 : 0.01568 : 0.017 : 0.01881 : 0.02291
  //6: 
  //7:
  //8:

  //with gaus bg
  //0: 4.6e5 : 8.754 : 4.5e5
  //1: 21.56 : 19 : 24.3
  //2: 3.4 : 3.13 : 3.5
  //3: 
  //4: 
  //5: 
  //6: 1.77e4, 4e4, 1.9e5
  //7: 8.9, 8.5, -19.55  
  //8: 9, 11.4, 14.89



  

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
