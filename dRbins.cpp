#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "RooCrystalBall.cxx"
#include "RooCrystalBall.h"
#include "MakeFitData.cpp"
#include "MakeFitMC.cpp"
#include "Efficiency.cpp"
#include <fstream>
#include <cmath>
#include <iomanip>
using std::setw;

void dRbins()
{
  int i;
  std::vector<std::pair<float,float>> YieldErrorDataPass, YieldErrorDataFail, YieldErrorMCPass, YieldErrorMCFail;
  float bins[]={0.0,0.15,0.2,0.25,0.3,0.35,0.5};  //this part is configurable, select bins for fitting
  bool FitInnerBarrel=1, FitOuterBarrel=0, FitEndcap=0; //this part is configurable, fit for barrel or endcap
  int NumberOfBins=sizeof(bins)/sizeof(bins[0])-1;
  TString InnerBarrelCut="fabs(ProbeEta)<0.8";
  TString OuterBarrelCut="fabs(ProbeEta)>0.8 && fabs(ProbeEta)<1.44";
  TString EndcapCut="fabs(ProbeEta)>1.57 && fabs(ProbeEta)<2.5";
  /*TFile *file_data = TFile::Open("TnPpairs_DATA.root");
  TTree *tree_data = (TTree*)file_data->Get("Events");
  tree_data->SetBranchStatus("*",0);
  tree_data->SetBranchStatus("Diele_mass",1);
  tree_data->SetBranchStatus("ele_pt",1);
  tree_data->SetBranchStatus("ProbePt",1);
  tree_data->SetBranchStatus("ProbeEta",1);
  tree_data->SetBranchStatus("ProbePass",1);
  gROOT->cd();*/

  TFile *file_MC = TFile::Open("TnPpairs_MC.root");
  TTree *tree_MC = (TTree*)file_MC->Get("Events");
  tree_MC->SetBranchStatus("*",0);
  tree_MC->SetBranchStatus("Diele_mass",1);
  tree_MC->SetBranchStatus("ele_pt",1);
  tree_MC->SetBranchStatus("ProbePt",1);
  tree_MC->SetBranchStatus("ProbeEta",1);
  tree_MC->SetBranchStatus("ProbePass",1);
  tree_MC->SetBranchStatus("dR",1);
  gROOT->cd();

  //Passing and failling probes (data and MC)
  //TTree* PassingProbesData[NumberOfBins];
  //TTree* FaillingProbesData[NumberOfBins];
  TTree* PassingProbesMC[NumberOfBins];
  TTree* FaillingProbesMC[NumberOfBins];

  if(FitInnerBarrel)
  for(i=0;i<NumberOfBins;i++)
  {
     //PassingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==1 && " + InnerBarrelCut,bins[i],bins[i+1]));
     //FaillingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==0 && " + InnerBarrelCut,bins[i],bins[i+1]));
     PassingProbesMC[i] = tree_MC->CopyTree(TString::Format("dR>%f && dR<%f && ProbePass==1 && ProbePt>5 && ProbePt<7 &&" + InnerBarrelCut,bins[i],bins[i+1]));
     FaillingProbesMC[i] = tree_MC->CopyTree(TString::Format("dR>%f && dR<%f && ProbePass==0 && ProbePt>5 && ProbePt<7 &&" + InnerBarrelCut,bins[i],bins[i+1]));
  }

  if(FitOuterBarrel)
  for(i=0;i<NumberOfBins;i++)
  {
     //PassingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==1 && " + OuterBarrelCut,bins[i],bins[i+1]));
     //FaillingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==0 && " + OuterBarrelCut,bins[i],bins[i+1]));
     PassingProbesMC[i] = tree_MC->CopyTree(TString::Format("dR>%f && dR<%f && ProbePass==1 && ProbePt>5 &&" + OuterBarrelCut,bins[i],bins[i+1]));
     FaillingProbesMC[i] = tree_MC->CopyTree(TString::Format("dR>%f && dR<%f && ProbePass==0 && ProbePt>5 &&" + OuterBarrelCut,bins[i],bins[i+1]));
  }

  if(FitEndcap)
  for(i=0;i<NumberOfBins;i++)
  {
     //PassingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==1 && " + EndcapCut,bins[i],bins[i+1]));
     //FaillingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==0 && " + EndcapCut,bins[i],bins[i+1]));
     PassingProbesMC[i] = tree_MC->CopyTree(TString::Format("dR>%f && dR<%f && ProbePass==1 && ProbePt>5 &&" + EndcapCut,bins[i],bins[i+1]));
     FaillingProbesMC[i] = tree_MC->CopyTree(TString::Format("dR>%f && dR<%f && ProbePass==0 && ProbePt>5 &&" + EndcapCut,bins[i],bins[i+1]));
  }

  //Fit for passing and failling probes in each bin
  for(i=0;i<NumberOfBins;i++)
  {
    //YieldErrorDataPass.push_back(MakeFitData(PassingProbesData[i], 1,  "plots/PassingProbesData_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf"));
    //YieldErrorDataFail.push_back(MakeFitData(FaillingProbesData[i], 1, "plots/FaillingProbesData_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf"));
    YieldErrorMCPass.push_back(MakeFitMC(PassingProbesMC[i], 1, "plots/PassingProbesMC_"+std::to_string(bins[i])+"_"+std::to_string(bins[i+1])+".pdf"));
    YieldErrorMCFail.push_back(MakeFitMC(FaillingProbesMC[i], 1, "plots/FaillingProbesMC_"+std::to_string(bins[i])+"_"+std::to_string(bins[i+1])+".pdf"));
  }

  TEfficiency* Eff_MC = Efficiency(YieldErrorMCPass, YieldErrorMCFail, bins, NumberOfBins);
  //TEfficiency* Eff_DATA = Efficiency(YieldErrorDataPass, YieldErrorDataFail, bins, NumberOfBins);

  //scale factors
  /*TH1F* SF = new TH1F("scale factor", "scale factor",NumberOfBins, bins);
  double SF_error_up[NumberOfBins], SF_error_dn[NumberOfBins];
  for(i=1;i<=NumberOfBins;i++)
  {
    SF->SetBinContent(i,Eff_DATA->GetEfficiency(i)*1.0/Eff_MC->GetEfficiency(i));
    SF_error_up[i]=sqrt(pow(Eff_MC->GetEfficiency(i)*Eff_DATA->GetEfficiencyErrorUp(i),2)+pow(Eff_DATA->GetEfficiency(i)*Eff_MC->GetEfficiencyErrorUp(i),2))/pow(Eff_MC->GetEfficiency(i),2);
    SF_error_dn[i]=sqrt(pow(Eff_MC->GetEfficiency(i)*Eff_DATA->GetEfficiencyErrorLow(i),2)+pow(Eff_DATA->GetEfficiency(i)*Eff_MC->GetEfficiencyErrorLow(i),2))/pow(Eff_MC->GetEfficiency(i),2);
    SF->SetBinError(i, SF_error_up[i]); //error up and down are practically equal so I put symmetric error equal to error up
  }
  SF->GetXaxis()->SetTitle("pT/GeV");
  SF->GetYaxis()->SetTitle("SF");
  TCanvas* SF_canvas = new TCanvas();
  gStyle->SetOptStat(0);
  SF->Draw("E");
  SF_canvas->SaveAs("plots/ScaleFactor.pdf");*/

  //efficiency plot
  //auto legend = new TLegend(0.9,0.9,1.0,1.0);
  auto mg  = new TMultiGraph();
  Eff_MC->SetLineColor(kRed);
  //Eff_DATA->SetLineColor(kBlue);
  TCanvas* oi = new TCanvas();
  oi->cd();
  Eff_MC->Draw();
  gPad->Update();
  auto graph_MC = Eff_MC->GetPaintedGraph();
  mg->Add(graph_MC);
  //Eff_DATA->Draw();
  //gPad->Update();
  //auto graph_data = Eff_DATA->GetPaintedGraph();
  //mg->Add(graph_data);
  mg->Draw("ap");
  mg->GetXaxis()->SetTitle("dR");
  mg->GetYaxis()->SetTitle("Efficiency");
  //legend->AddEntry(graph_data,"data");
  //legend->AddEntry(graph_MC,"MC");
  //legend->Draw();
  oi->SaveAs("plots/Efficiency.pdf");

  ofstream MCeff;
  MCeff.open ("eff_MC_dR.txt");
  MCeff <<std::left<<std::setw(9)<< "dR"<<std::setw(13)<<"eff"<<std::setw(13)<<"eff_err_up"<<"eff_error_dn"<<std::endl;
  for(i=1;i<=NumberOfBins;i++)
  MCeff<<std::setw(2)<<bins[i-1]<<"-"<<std::setw(6)<<bins[i]<<std::setw(10)<<Eff_MC->GetEfficiency(i)<<"   "<<std::setw(13)<<Eff_MC->GetEfficiencyErrorUp(i)<<std::setw(10)<<Eff_MC->GetEfficiencyErrorLow(i)<<std::endl;


  //Write into txt file
  /*ofstream SFfile;
  SFfile.open ("SF.txt");
  SFfile <<std::left<<std::setw(9)<< "pT"<<std::setw(13)<<"SF"<<std::setw(13)<<"SF_err_up"<<"SF_err_dn"<<std::endl;
  for(i=1;i<=NumberOfBins;i++)
  SFfile<<std::setw(2)<<bins[i-1]<<"-"<<std::setw(6)<<bins[i]<<std::setw(10)<<SF->GetBinContent(i)<<"   "<<std::setw(13)<<SF_error_up[i]<<std::setw(10)<<SF_error_dn[i]<<std::endl;

  ofstream MCeff;
  MCeff.open ("eff_MC.txt");
  MCeff <<std::left<<std::setw(9)<< "pT"<<std::setw(13)<<"eff"<<std::setw(13)<<"eff_err_up"<<"eff_error_dn"<<std::endl;
  for(i=1;i<=NumberOfBins;i++)
  MCeff<<std::setw(2)<<bins[i-1]<<"-"<<std::setw(6)<<bins[i]<<std::setw(10)<<Eff_MC->GetEfficiency(i)<<"   "<<std::setw(13)<<Eff_MC->GetEfficiencyErrorUp(i)<<std::setw(10)<<Eff_MC->GetEfficiencyErrorLow(i)<<std::endl;

  ofstream DATAeff;
  DATAeff.open ("eff_DATA.txt");
  DATAeff <<std::left<<std::setw(9)<< "pT"<<std::setw(13)<<"eff"<<std::setw(13)<<"eff_err_up"<<"eff_error_dn"<<std::endl;
  for(i=1;i<=NumberOfBins;i++)
  DATAeff<<std::setw(2)<<bins[i-1]<<"-"<<std::setw(6)<<bins[i]<<std::setw(10)<<Eff_DATA->GetEfficiency(i)<<"   "<<std::setw(13)<<Eff_DATA->GetEfficiencyErrorUp(i)<<std::setw(10)<<Eff_DATA->GetEfficiencyErrorLow(i)<<std::endl;
  */
}
