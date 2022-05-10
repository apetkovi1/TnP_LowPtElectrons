#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "RooCrystalBall.cxx"
#include "RooCrystalBall.h"
#include "MakeFitData.cpp"
#include "MakeFitMC.cpp"
#include "Efficiency.cpp"

void TnPAnalyzer()
{
  int i;
  std::vector<std::pair<float,float>> YieldErrorDataPass, YieldErrorDataFail, YieldErrorMCPass, YieldErrorMCFail;
  float bins[]={7,9,11,13,15};  //this part is configurable, select bins for fitting
  bool FitBarrel=1, FitEndcap=0; //this part is configurable, fit for barrel or endcap
  int NumberOfBins=sizeof(bins)/sizeof(bins[0])-1;
  TString BarrelCut="ProbeEta>-1.45 && ProbeEta<1.45";
  TString EndcapCut="((ProbeEta>-2.5 && ProbeEta<-1.45) || (ProbeEta>1.45 && ProbeEta<2.5))";
  TFile *file_data = TFile::Open("TnPpairs_DATA.root");
  TTree *tree_data = (TTree*)file_data->Get("Events");
  tree_data->SetBranchStatus("*",0);
  tree_data->SetBranchStatus("Diele_mass",1);
  tree_data->SetBranchStatus("ele_pt",1);
  tree_data->SetBranchStatus("mvaEleID_Fall17_noIso_V2_wp90",1);
  tree_data->SetBranchStatus("ProbePt",1);
  tree_data->SetBranchStatus("ProbeEta",1);
  tree_data->SetBranchStatus("ProbePass",1);
  gROOT->cd();

  TFile *file_MC = TFile::Open("TnPpairs_MC.root");
  TTree *tree_MC = (TTree*)file_MC->Get("Events");
  tree_MC->SetBranchStatus("*",0);
  tree_MC->SetBranchStatus("Diele_mass",1);
  tree_MC->SetBranchStatus("ele_pt",1);
  tree_MC->SetBranchStatus("mvaEleID_Fall17_noIso_V2_wp90",1);
  tree_MC->SetBranchStatus("ProbePt",1);
  tree_MC->SetBranchStatus("ProbeEta",1);
  tree_MC->SetBranchStatus("ProbePass",1);
  gROOT->cd();

  //Passing and failling probes (data and MC)
  TTree* PassingProbesData[NumberOfBins];
  TTree* FaillingProbesData[NumberOfBins];
  TTree* PassingProbesMC[NumberOfBins];
  TTree* FaillingProbesMC[NumberOfBins];

  if(FitBarrel)
  for(i=0;i<NumberOfBins;i++)
  {
     PassingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==1 && " + BarrelCut,bins[i],bins[i+1]));
     FaillingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==0 && " + BarrelCut,bins[i],bins[i+1]));
     PassingProbesMC[i] = tree_MC->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==1 && " + BarrelCut,bins[i],bins[i+1]));
     FaillingProbesMC[i] = tree_MC->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==0 && " + BarrelCut,bins[i],bins[i+1]));
  }

  if(FitEndcap)
  {
     PassingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==1 && " + EndcapCut,bins[i],bins[i+1]));
     FaillingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==0 && " + EndcapCut,bins[i],bins[i+1]));
     PassingProbesMC[i] = tree_MC->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==1 && " + EndcapCut,bins[i],bins[i+1]));
     FaillingProbesMC[i] = tree_MC->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==0 && " + EndcapCut,bins[i],bins[i+1]));
  }

  //Fit for passing and failling probes in each bin
  for(i=0;i<NumberOfBins;i++)
  {
    YieldErrorDataPass.push_back(MakeFitData(PassingProbesData[i], 1,  "plots/PassingProbesData_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf"));
    YieldErrorDataFail.push_back(MakeFitData(FaillingProbesData[i], 1, "plots/FaillingProbesData_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf"));
    YieldErrorMCPass.push_back(MakeFitMC(PassingProbesMC[i], 1, "plots/PassingProbesMC_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf"));
    YieldErrorMCFail.push_back(MakeFitMC(FaillingProbesMC[i], 1, "plots/FaillingProbesMC_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf"));
  }

  Efficiency(YieldErrorDataPass,YieldErrorDataFail,YieldErrorMCPass, YieldErrorMCFail, bins, NumberOfBins);
}
