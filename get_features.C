//////////////////////////////////////////////////////////
//////// Get Features of Historam
////////    Finds local minima, maxima, dydx=0 and d2ydx2=0 point
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

#include "Fit/FitResult.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "TRandom3.h"
#include "TDatime.h"

#include "Math/MinimizerOptions.h"

using namespace std;

bool b_draw = true;//false; // draw output histograms?
  
void get_features()
{
  TH1D * h_minima   = new TH1D("h_minima","h_minima",50,0,50);
  TH1D * h_maxima   = new TH1D("h_maxima","h_maxima",50,0,50);
  TH1D * h_dxdy0    = new TH1D("h_dxdy0","h_dxdy0",50,0,50);
  TH1D * h_d2ydx2_0 = new TH1D("h_d2ydx2_0","h_d2ydx2_0",50,0,50);

  // Write Features to txt file for R
  ofstream ofs_out;
  //ofs_out.open("out_all.txt");
  ofs_out.open("demo.txt");
  
  // find all minima and maxima
  TFile * tf_in = new TFile("training_data.root","READ");

  // Number of historgrams to read
  int i_N_histos = 100000;
  if(b_draw) i_N_histos = 2;
  if(i_N_histos == 0) return;

  // set fit function
  const char *la0  = "[3]*TMath::Landau(x,[1],[2],0)";
  const char *la1  = "(1-[3])*TMath::Landau(x,2*[1]+1.4*[2],2.0*[2],0)";
  const char *c_bg = "[4]*exp(-[5]*x)+[6]*exp(-0.5*((x-[7])/[8])**2)";

  TFormula * tform_s_bg = new TFormula("tform_s_bg", Form("[0]*(%s+%s)+%s",la0,la1,c_bg));

  bool b_draw_derivatives_on_main_plot = true;//false;

  for(int ijk = 0; ijk < i_N_histos; ijk++)
    {
      if(b_draw&&(ijk==0)) ijk+=2;
      // Load Histogram from data file
      TH1D * h_energy = (TH1D*) tf_in -> Get(Form("h_energy_%d",ijk));
      if(!h_energy) return;
      h_energy -> SetXTitle("Charged Deposited in Sensor");
      h_energy -> SetYTitle("Counts");


      //Get the Generating Function, for the training response
      TF1 * f1 = (TF1*) tf_in -> Get(Form("f1_%d",ijk));

      // Get Response MIP MPV
      double d_training_response = f1 -> GetParameter(1);
      
      // Draw Histogram
      TCanvas * TC_energy = new TCanvas(Form("TC_energy_%d",ijk),"TC_energy");
      TC_energy -> cd();
      gPad -> SetLogy(); gPad->SetGridx(); gPad->SetGridy();  h_energy->Draw("HIST");// f1 -> Draw("same");

      // make a legend
      TLegend * TL_legend = new TLegend(0.55,0.525,0.99,0.935);
      //            TLegend * TL_legend = new TLegend(0.55,0.525,0.85,0.935);
      //TL_legend -> AddEntry(h_energy,"SensorEnergy Deposition");


      //////////////////////////////////////////////////////////
      /////// Local Minima and Maxima
      //////////////////////////////////////////////////////////      
      double  d_adc_start = 80;
      double  d_adc_end   = 0;
      int iBin_adc_start = h_energy -> FindFixBin(d_adc_start);
      int iBin_adc_end   = h_energy -> FindFixBin(d_adc_end);      
      int index = 0;
      double d_tmp_min = -9999.0;
      double d_tmp_max = -9999.0;
      bool b_pos_slope = 0;

      //record the maximum closest to 0
      int i_max0 = -9999;
      double d_max0_dist = -9999.0;

      vector<int> v_iBin_minima;
      vector<int> v_iBin_maxima;
      vector<double> v_pos_maxima;
      vector<double> v_pos_minima;
      for(int iBin = 1; iBin < iBin_adc_start; iBin ++)
	{
	  double d_center  = h_energy -> GetBinCenter(iBin);
	  double d_content = h_energy -> GetBinContent(iBin);

	  bool b_min = true;
	  bool b_max = true;

	  double a5_pre_content[5] = {0.,0.,0.,0.,0.};
	  double a5_post_content[5] = {0.,0.,0.,0.,0.};

	  int i_num_check = 5;
	  for(int jkl = 0; jkl < i_num_check;jkl++)
	    {	      
	      if((iBin-1-jkl) > 0) a5_pre_content[jkl]  = h_energy -> GetBinContent(iBin-1-jkl);
	      else { b_min = false; b_max = false;}
	      if((iBin+1+jkl) < h_energy ->GetNbinsX()) a5_post_content[jkl] = h_energy -> GetBinContent(iBin+1+jkl);
	      else { b_min = false; b_max = false;}

	      if(a5_pre_content[jkl] < d_content) b_min = false;
	      if(a5_pre_content[jkl] > d_content) b_max = false;
		
	      if(a5_post_content[jkl] < d_content) b_min = false;
	      if(a5_post_content[jkl] > d_content) b_max = false;

	      if(jkl==0) continue;
	      if(a5_pre_content[jkl]  == d_content) b_max = false;
	      if(a5_post_content[jkl] == d_content) b_min = false;   
		
	    }//for(int jkl = 0; jkl < 5;jkl++)
	      
	  double d_pre_center  = h_energy -> GetBinCenter(iBin-1);
	  double d_pre_content = h_energy -> GetBinContent(iBin-1);
	  int iBin_pre = iBin-1;

	  double d_next_center  = h_energy -> GetBinCenter(iBin+1);
	  double d_next_content = h_energy -> GetBinContent(iBin+1);
	  int iBin_next = iBin+1;
    
	  if(iBin == iBin_adc_start) continue;
	    
	  if(b_max)
	    {
	      v_iBin_maxima.push_back(iBin);
	      v_pos_maxima.push_back(d_center);
	      if( fabs(d_center-0) < fabs(d_max0_dist) )
		{
		  i_max0 = iBin;
		  d_max0_dist = fabs(d_center-0);
		}//if( fabs(d_center-0) < fabs(d_max0_dist) )
	    }//if((d_pre_content < d_content) && (d_next_content < d_content) )
	  if(b_min)
	    { 
	      v_iBin_minima.push_back(iBin);
	      v_pos_minima.push_back(d_center);
	      //cout<<"          Minima: "<<d_center<<endl; 
	    } //if((d_pre_content > d_content) && (d_next_content > d_content) ) 
	}//for(int iBin = iBin_adc_start; iBin > iBin_adc_end; iBin --)
      
      double d_mean = h_energy -> GetBinCenter(i_max0);
      double d_max  = h_energy -> GetBinContent(i_max0);

      int i_min0 = -9999;
      double d_min0_dist = -9999;   
      vector<TArrow*> a_ptr_arrows;
      for(int k = 0; k< v_iBin_minima.size(); k++)
	{
	  double d_center  = h_energy -> GetBinCenter(v_iBin_minima[k]);
	  double d_content = h_energy -> GetBinContent(v_iBin_minima[k]);
	  // TLine * TL_minima = new TLine(d_center,0.0,d_center,h_energy->GetMaximum());
	  // TL_minima -> SetLineColor(kRed);
	  // TL_minima -> Draw();

	  TArrow * TL_minima = new TArrow(d_center,d_content-(d_content*0.1),d_center,d_content,0.015,"|>");
	  TL_minima -> SetLineColor(kGreen);
	  TL_minima -> SetFillColor(kGreen);

	  TL_minima -> Draw();

	  if(k==0) TL_legend -> AddEntry(TL_minima ,"Minima");
	  
	  a_ptr_arrows.push_back(TL_minima);
	  //find the minimum closest to the candidate mean
	  if( fabs(d_center-d_mean) < fabs(d_min0_dist) )
	    {
	      i_min0 = k;
	      //cout<<" -"<<d_center<<" "<<d_mean<<endl;
	      d_min0_dist = fabs(d_center-d_mean);
	    }//if( fabs(d_center-0) < fabs(d_min0_dist) )	  
	}//for(int k = 0; k< v_iBin_minima.size(); k++)

      double d_min = h_energy -> GetBinContent(i_min0);
      
      for(int k = 0; k< v_iBin_maxima.size(); k++)
	{
	  double d_center  = h_energy -> GetBinCenter(v_iBin_maxima[k]);
	  double d_content = h_energy -> GetBinContent(v_iBin_maxima[k]);
	  // TLine * TL_maxima = new TLine(d_center,0.0,d_center,h_energy->GetMaximum());
	  // TL_maxima -> SetLineColor(kBlue);
	  // TL_maxima -> Draw();

	  TArrow * TL_maxima = new TArrow(d_center,d_content-(d_content*0.1),d_center,d_content,0.015,"|>");
	  TL_maxima -> SetLineColor(kBlue);
	  TL_maxima -> SetFillColor(kBlue);
	  TL_maxima -> Draw();
	  a_ptr_arrows.push_back(TL_maxima);
	  if(k==0) TL_legend -> AddEntry(TL_maxima ,"Maxima");
	}//for(int k = 0; k< v_iBin_maxima.size(); k++)

      TArrow * TL_mean = new TArrow(d_mean,d_max,d_mean,d_max+(d_max*0.1),0.015,"|>");
      TL_mean -> SetLineColor(kBlue);
      TL_mean -> SetFillColor(kBlue);
      TL_mean -> Draw();
      a_ptr_arrows.push_back(TL_mean);

      //////////////////////////////////////////////////////////
      /// END Local Minima and Maxima
      //////////////////////////////////////////////////////////

      ////////////////////////////////////////////////
      /// /// local dervatives
      ////////////////////////////////////////////////

      // TH1D * h_dydx = (TH1D*) h_energy -> Clone(Form("h_dydx_%d",ijk));
      // h_dydx -> Reset();
      TH1D * h_dydx = new TH1D(Form("h_dydx_%d",ijk),"dydx",150,0,150);      
      for(int iBin = iBin_adc_start; iBin > iBin_adc_end; iBin --)
	{
	  if((iBin <2)||(iBin > (h_energy->GetNbinsX()+1)) ) continue;
	  double dx = h_energy -> GetBinWidth(1);
	  double d_center  = h_energy -> GetBinCenter(iBin);
	  double d_content = h_energy -> GetBinContent(iBin);

	  double d_pre_content  = h_energy -> GetBinContent(iBin-11);
	  double d_next_content  = h_energy -> GetBinContent(iBin+1);

	  //	    double dydx = (d_next_content-d_pre_content) / (2*dx);
	  //	    double dydx = (d_content-d_pre_content) / (dx);
	  double dydx = (d_next_content-d_content) / (dx);
	  // h_dydx -> Fill(dydx);
	  h_dydx -> SetBinContent(iBin,dydx);
	    
	}//for(int iBin = iBin_adc_start; iBin > iBin_adc_end; iBin --)
      TCanvas * TC_dydx = new TCanvas(Form("TC_dydx_%d",ijk),"TC_dydx");
      TC_dydx -> cd();
      gPad->SetGridx(); gPad->SetGridy();
      h_dydx -> Smooth();
      h_dydx -> Draw("HIST");

      // find dydx = 0
      vector<int> v_iBin_dxdy0;
      for(int iBin = 1; iBin < (h_dydx -> GetNbinsX()-1); iBin++)
	{
	  bool b_dydx0 = false;
	  double d_center  = h_dydx -> GetBinCenter(iBin);
	  double d_content = h_dydx -> GetBinContent(iBin);

	  double d_pre_content  = h_dydx -> GetBinContent(iBin-1);
	  double d_next_content  = h_dydx -> GetBinContent(iBin+1);  

	  if((d_content < 0) && (d_pre_content > 0)) b_dydx0 = true;
	  if((d_content > 0) && (d_pre_content < 0)) b_dydx0 = true;

	  // if((d_content < 0) && (d_next_content > 0)) b_dydx0 = true;
	  // if((d_content > 0) && (d_next_content < 0)) b_dydx0 = true;

	  // if((d_pre_content < 0) && (d_content > 0)) b_dydx0 = true;
	  // if((d_pre_content > 0) && (d_content < 0)) b_dydx0 = true;

	  if(!b_dydx0) continue;
	  // cout<<" dxdy0: "<<d_pre_content
	  //     <<" "<<d_next_content
	  //     <<endl;

	  if(b_dydx0) v_iBin_dxdy0.push_back(d_center);
	}//for(int iBin = 1; iBin < (h_dydx -> GetNbinsX()-1); iBin++)

      for(int jkl = 0; jkl < v_iBin_dxdy0.size(); jkl++)
	{
	  double d_center    = h_dydx -> GetBinCenter(v_iBin_dxdy0[jkl]);
	  double d_content   = h_dydx -> GetBinContent(v_iBin_dxdy0[jkl]);

	  int iBin = h_energy -> FindFixBin(v_iBin_dxdy0[jkl]);								     
	  double d_Econtent = h_energy -> GetBinContent(iBin);
	  
	  TArrow * TL_dxdy0 = new TArrow(d_center,d_content-(d_content*0.1),d_center,d_content,0.015,"|>");
	  TL_dxdy0 -> SetLineColor(kRed);
	  TL_dxdy0 -> SetFillColor(kRed);	  

	  TC_dydx  -> cd();
	  TL_dxdy0 -> Draw();

	  //	  if(jkl==0) TL_legend -> AddEntry(TL_dxdy0 ,"#frac{dy}{dx}=0");
	  if(jkl==0) TL_legend -> AddEntry(TL_dxdy0 ,"dy/dx=0");

	  if(b_draw_derivatives_on_main_plot)
	    {
	      TC_energy  -> cd();
	      //TL_dxdy0 -> DrawArrow(d_center+1,d_Econtent-(d_Econtent*0.1),d_center+1,d_Econtent);
	      TL_dxdy0 -> DrawArrow(d_center-1,d_Econtent-(d_Econtent*0.1),d_center-1,d_Econtent);
	    }//if(b_draw_derivatives_on_main_plot)

	  a_ptr_arrows.push_back(TL_dxdy0);
	}//for(int jkl = 0; jkl < v_iBin_dxdy0.size(); jkl++)

      ////////////////////////////////////////////////
      /// END local dervatives
      ////////////////////////////////////////////////

      ////////////////////////////////////////////////
      /// /// local 2nd dervatives
      ////////////////////////////////////////////////
      TH1D * h_d2ydx2 = new TH1D(Form("h_d2ydx2_%d",ijk),"d2ydx2",150,0,150);      
      for(int iBin = iBin_adc_start; iBin > iBin_adc_end; iBin --)
	{
	  if((iBin <2)||(iBin > (h_dydx->GetNbinsX()+1)) ) continue;
	  double dx = h_dydx -> GetBinWidth(1);
	  double d_center  = h_dydx -> GetBinCenter(iBin);
	  double d_content = h_dydx -> GetBinContent(iBin);

	  double d_pre_content  = h_dydx -> GetBinContent(iBin-11);
	  double d_next_content  = h_dydx -> GetBinContent(iBin+1);

	  //	    double d2ydx2 = (d_next_content-d_pre_content) / (2*dx);
	  //	    double d2ydx2 = (d_content-d_pre_content) / (dx);
	  double d2ydx2 = (d_next_content-d_content) / (dx);
	  // h_d2ydx2 -> Fill(d2ydx2);
	  h_d2ydx2 -> SetBinContent(iBin,d2ydx2);	    
	}//for(int iBin = iBin_adc_start; iBin > iBin_adc_end; iBin --)

      TCanvas * TC_d2ydx2 = new TCanvas(Form("TC_d2ydx2_%d",ijk),"TC_d2ydx2");
      TC_d2ydx2 -> cd();
      gPad->SetGridx(); gPad->SetGridy();
      h_d2ydx2 -> Smooth();
      h_d2ydx2 -> Draw("HIST");

      // find d2ydx2 = 0
      vector<int> v_iBin_d2ydx2_0;
      for(int iBin = 1; iBin < (h_dydx -> GetNbinsX()-1); iBin++)
	{
	  bool b_d2ydx2_0 = false;
	  double d_center  = h_d2ydx2 -> GetBinCenter(iBin);
	  double d_content = h_d2ydx2 -> GetBinContent(iBin);

	  double d_pre_content  = h_d2ydx2 -> GetBinContent(iBin-1);
	  double d_next_content  = h_d2ydx2 -> GetBinContent(iBin+1);  

	  if((d_content < 0) && (d_pre_content > 0)) b_d2ydx2_0 = true;
	  if((d_content > 0) && (d_pre_content < 0)) b_d2ydx2_0 = true;

	  // if((d_content < 0) && (d_next_content > 0)) b_d2ydx2_0 = true;
	  // if((d_content > 0) && (d_next_content < 0)) b_d2ydx2_0 = true;

	  // if((d_pre_content < 0) && (d_content > 0)) b_d2ydx2_0 = true;
	  // if((d_pre_content > 0) && (d_content < 0)) b_d2ydx2_0 = true;

	  if(!b_d2ydx2_0) continue;
	  // cout<<" d2ydx2_0: "<<d_pre_content
	  //     <<" "<<d_next_content
	  //     <<endl;

	  if(b_d2ydx2_0) v_iBin_d2ydx2_0.push_back(d_center);
	}//for(int iBin = 1; iBin < (h_d2ydx2 -> GetNbinsX()-1); iBin++)

      for(int jkl = 0; jkl < v_iBin_d2ydx2_0.size(); jkl++)
	{
	  double d_center    = h_d2ydx2 -> GetBinCenter(v_iBin_d2ydx2_0[jkl]);
	  double d_content   = h_d2ydx2 -> GetBinContent(v_iBin_d2ydx2_0[jkl]);

	  int iBin = h_energy -> FindFixBin(v_iBin_d2ydx2_0[jkl]);								     
	  double d_Econtent = h_energy -> GetBinContent(iBin);
	  
	  TArrow * TL_d2ydx2_0 = new TArrow(d_center,d_content-(d_content*0.1),d_center,d_content,0.015,"|>");
	  TL_d2ydx2_0 -> SetLineColor(kBlack);
	  TL_d2ydx2_0 -> SetFillColor(kBlack);

	  TC_d2ydx2  -> cd();
	  TL_d2ydx2_0 -> Draw();

	  //if(jkl==0) TL_legend -> AddEntry(TL_d2ydx2_0 ,"#frac{d^{2}y}{dx^{2}}=0");
	  if(jkl==0) TL_legend -> AddEntry(TL_d2ydx2_0 ,"d^{2}y/dx^{2}=0");

	  if(b_draw_derivatives_on_main_plot)
	    {
	      TC_energy  -> cd();
	      TL_d2ydx2_0 -> DrawArrow(d_center-1,d_Econtent-(d_Econtent*0.1),d_center-1,d_Econtent);
	    }//if(b_draw_derivatives_on_main_plot)
	  a_ptr_arrows.push_back(TL_d2ydx2_0);
	}//for(int jkl = 0; jkl < v_iBin_d2ydx2_0.size(); jkl++)
      ////////////////////////////////////////////////
      /// END local 2nd dervatives
      ////////////////////////////////////////////////

      h_minima   ->Fill(v_iBin_minima.size());
      h_maxima   ->Fill(v_iBin_maxima.size());
      h_dxdy0    ->Fill(v_iBin_dxdy0.size());
      h_d2ydx2_0 ->Fill(v_iBin_d2ydx2_0.size());

      // cout<<" Num Minima: "<<v_iBin_minima.size()   <<endl;
      // cout<<" Num Maxima: "<<v_iBin_maxima.size()   <<endl;
      // cout<<" Num dydx:   "<<v_iBin_dxdy0.size()    <<endl;
      // cout<<" Num d2ydx2: "<<v_iBin_d2ydx2_0.size() <<endl;


      TC_dydx   -> Close();
      TC_d2ydx2 -> Close();

      //    Response: location of mpv
      // Vars:
      // xpos:    max 0, 1, 2, 3
      //          dydx 0, 1, 2, 3
      //          dy2dx 0, 1, 2, 3
      // dist: max_i from dydx
      //       max_i from d2ydx2
      ///  OR ORDER?
      // fit each max,dydx,d2ydx2 to landau and use chi2 as variable

      //output vectors
      vector<double> v_out_min;	  
      vector<double> v_out_max;	  
      vector<double> v_out_dydx;  
      vector<double> v_out_d2ydx2;

      vector<vector<double>> v2_dist_max_dydx;
      vector<vector<double>> v2_dist_max_d2ydx2;
      //      for(int jkl = 0; jkl<20; jkl++)
      for(int jkl = 0; jkl<15; jkl++)
	{
	  if(jkl < 5)
	    {
	      if( jkl < v_iBin_minima.size() ) v_out_min.push_back(v_iBin_minima[jkl]);
	      else v_out_min.push_back(-9999.0);

	      vector<double>  v_dist_max_dydx;
	      vector<double>  v_dist_max_d2ydx2;
	      if( jkl < v_iBin_maxima.size() )
		{
		  v_out_max.push_back(v_iBin_maxima[jkl]);
	      
		  for(int klm = 0; klm < 14; klm++)
		    {
		      if(klm<7)
			{
			  if(klm < v_iBin_dxdy0.size() ) v_dist_max_dydx.push_back(v_iBin_maxima[jkl]-v_iBin_dxdy0[klm]);
			  else v_dist_max_dydx.push_back(-9999.0);
			}//if(klm<7)
		      
		      if(klm < v_iBin_d2ydx2_0.size() ) v_dist_max_d2ydx2.push_back(v_iBin_maxima[jkl]-v_iBin_d2ydx2_0[klm]);
		      else v_dist_max_d2ydx2.push_back(-9999.0);
		    }//for(int klm = 0; klm < 14; klm++)
		}//if( jkl < v_iBin_maxima.size() )
	      else
		{
		  v_out_max.push_back(-9999.0);
		  for(int klm = 0; klm < 14; klm++)
		    {
		      if(klm<7)			
			{
			  v_dist_max_dydx.push_back(-9999.0);
			}//if(klm<7)
		      v_dist_max_d2ydx2.push_back(-9999.0);
		    }//for(int klm = 0; klm < 14; klm++)
		}//else( jkl < v_iBin_maxima.size() )

	      v2_dist_max_dydx.push_back(v_dist_max_dydx);
	      v2_dist_max_d2ydx2.push_back(v_dist_max_d2ydx2);	      
	    }//if(jkl < 5)
	  if(jkl < 7)
	    {
	      if( jkl < v_iBin_dxdy0.size() ) v_out_dydx.push_back(v_iBin_dxdy0[jkl]);
	      else v_out_dydx.push_back(-9999.0);
	    }//if(jkl < 7)
	  if( jkl < v_iBin_d2ydx2_0.size() ) v_out_d2ydx2.push_back(v_iBin_d2ydx2_0[jkl]);
	  else v_out_d2ydx2.push_back(-9999.0);

	}//for(int jkl = 0; jkl<20;jkl++)      

      //write to output txt file
      if( ijk == 0)
	{
	  for(int jkl = 0; jkl<v_out_min.size(); jkl++) ofs_out<<"Min"<<jkl<<" ";
	  for(int jkl = 0; jkl<v_out_max.size(); jkl++) ofs_out<<"Max"<<jkl<<" ";
	  for(int jkl = 0; jkl<v_out_dydx.size(); jkl++) ofs_out<<"dydx"<<jkl<<" ";
	  for(int jkl = 0; jkl<v_out_d2ydx2.size(); jkl++) ofs_out<<"d2ydx2_"<<jkl<<" ";

	  for(int jkl = 0; jkl<v_out_max.size(); jkl++)
	    {
	      for(int klm = 0; klm < v2_dist_max_dydx[jkl].size(); klm++)
		{	      
		  ofs_out<< "dist_"<<jkl<<"_"<<klm<<" ";
		}//for(int klm = 0; klm < v2_dist_max_dydx[jkl].size(); klm++)

	      for(int klm = 0; klm < v2_dist_max_d2ydx2[jkl].size(); klm++)
		{	      
		  ofs_out<< "dist2_"<<jkl<<"_"<<klm<<" ";
		}//for(int klm = 0; klm < v2_dist_max_d2ydx2[jkl].size(); klm++)
	    }//for(int jkl = 0; jkl<v_out_max.size(); jkl++)
	  ofs_out<<" MPV"<<endl;		
	}//if( ijk == 0)
      // data out txt
      for(int jkl = 0; jkl<v_out_min.size(); jkl++)    ofs_out<<v_out_min[jkl]<<" ";
      for(int jkl = 0; jkl<v_out_max.size(); jkl++)    ofs_out<<v_out_max[jkl]<<" ";
      for(int jkl = 0; jkl<v_out_dydx.size(); jkl++)   ofs_out<<v_out_dydx[jkl]<<" ";
      for(int jkl = 0; jkl<v_out_d2ydx2.size(); jkl++) ofs_out<<v_out_d2ydx2[jkl]<<" ";
      for(int jkl = 0; jkl<v_out_max.size(); jkl++)
	{
	  for(int klm = 0; klm < v2_dist_max_dydx[jkl].size(); klm++)
	    {
	      ofs_out<< v2_dist_max_dydx[jkl][klm] <<" ";
	    }//for(int klm = 0; klm < v2_dist_max_dydx[jkl].size(); klm++)	 
	  for(int klm = 0; klm < v2_dist_max_d2ydx2[jkl].size(); klm++)
	    {
	      ofs_out<< v2_dist_max_d2ydx2[jkl][klm] <<" ";
	    }//for(int klm = 0; klm < v2_dist_max_d2ydx2[jkl].size(); klm++)	  
	}//for(int jkl = 0; jkl<v_out_max.size(); jkl++)
      //record training response:
      ofs_out<<d_training_response;
      //for now record respones as intergers for easier classification
      //      ofs_out<< TMath::Floor(d_training_response);
      ofs_out<<endl;

      //delete
      if(!b_draw)
	{
	  delete h_energy;
	  delete f1;
	  delete TC_energy;
	  delete h_dydx;
	  delete TC_dydx;
	  delete h_d2ydx2;
	  delete TC_d2ydx2;

	  for(int jkl = 0; jkl < a_ptr_arrows.size(); jkl++)
	    {
	      delete  a_ptr_arrows[jkl];
	    }//for(int jkl = 0; jkl < a_ptr_arrows.size(); jkl++)
	}//if(!b_draw)
      else
	{
	        TL_legend -> Draw();
	}// else (!b_draw)
      //cout<<" ================================ "<<endl;
    }//for(int ijk = 0; ijk < i_N_histos; ijk++)

  // new TCanvas();
  // h_minima     -> Draw();
  // new TCanvas();
  // h_maxima     -> Draw();
  // new TCanvas();
  // h_dxdy0      -> Draw();
  // new TCanvas();
  // h_d2ydx2_0   -> Draw();

  return;
}//void get_features()

//is only read by compiler
#ifndef __CINT__

int main(int argc, char * argv[]){
  gSystem->Load("libTree");
  TApplication App("Analysis", &argc, argv);

  if(!b_draw) gROOT->SetBatch(kTRUE);
  get_features();
  
  //App.Terminate();
  cout << "Main: program exiting." << endl;
  pthread_exit(NULL);

  return 0;
}
#endif
