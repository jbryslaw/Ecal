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
#include "TArrow.h"

#include "TMultiLayerPerceptron.h"
#include "TMLPAnalyzer.h"


#include "Fit/FitResult.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "TRandom3.h"
#include "TDatime.h"

#include "Math/MinimizerOptions.h"

using namespace std;

void find_mpv_NN()
{
  TFile * tf_in = new TFile("out_all_bins_ttree.root","READ");
  TTree * t_traintest = (TTree*) tf_in -> Get("t_traintest");
  
  //set hidden layers
  vector<int> v_hidden_ly;
  v_hidden_ly.push_back(10);
  v_hidden_ly.push_back(8);

  //define objects to be used in neural net
  TString ts_percept_def;
  //  for(int ijk = 2; ijk < 81; ijk++)
  // just use bins 9 to 40, other bins are mostly irrelevant
  for(int ijk = 9; ijk < 40; ijk++)
    {
      if(ijk == 9) ts_percept_def = Form("BinContent%d",ijk);
      else ts_percept_def += Form(",BinContent%d",ijk);
    }//for(int ijk = 2; ijk < 40; ijk++)
  //set hidden layers
  for(int ijk = 0; ijk < v_hidden_ly.size(); ijk++)
    ts_percept_def += Form(":%d",v_hidden_ly[ijk]);
  // set output
  //  ts_percept_def += ":MPV";
  ts_percept_def += ":MPV";
  
  TMultiLayerPerceptron* TMLP_ecal_percept
    = new TMultiLayerPerceptron(ts_percept_def,t_traintest,"Entry$%2","(Entry$%2)==0");

  // Train (epoch == 1back and forth trough all training 
  TMLP_ecal_percept->Train(150,"graph update=10");

  // Analyze
  TMLPAnalyzer* TMLP_ecal_ana=new TMLPAnalyzer(TMLP_ecal_percept);
  TMLP_ecal_ana->GatherInformations();
  TMLP_ecal_ana->CheckNetwork();
  TMLP_ecal_ana->DrawDInputs();

  // histogram difference between test values and predictionvalues
  TH1D * h_diff_pred_true = new TH1D("h_diff_pred_true","Predicted MPV - True MPV",5,0,5);    
  TMLP_ecal_ana->GetIOTree()->Draw("fabs(Out.Out0-True.True0)>>h_diff_pred_true");

  //Calculate fraction of successful predictions
  double d_total_integral = h_diff_pred_true -> Integral(0,-1);  
  double d_less1_integral = h_diff_pred_true -> Integral(1,1);  
  double d_entries_less1  = h_diff_pred_true -> GetBinContent(1);
  double d_frac_less1     = (d_total_integral != 0.0) ? d_less1_integral/d_total_integral : -9999.0;
  cout<<" Fraction of Successful Predictions: "<<d_frac_less1
    //      <<" : "<<d_less1_integral<<" : "<<d_frac_less1<<" : "<<d_entries_less1<<" : "<<d_total_integral
      <<endl;

  return;
}//void find_mpv_NN()

//is only read by compiler
#ifndef __CINT__

int main(int argc, char * argv[]){
  gSystem->Load("libTree");
  TApplication App("Analysis", &argc, argv);

  find_mpv_NN();
  // gROOT->SetBatch(kTRUE);
  
  //App.Terminate();
  cout << "Main: program exiting." << endl;
  pthread_exit(NULL);

  return 0;
}
#endif
