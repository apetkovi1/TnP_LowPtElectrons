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
  float bins[]={5,7,10,15};  //this part is configurable, select bins for fitting
  bool FitInnerBarrel=1, FitOuterBarrel=0, FitEndcap=0; //this part is configurable, fit for barrel or endcap
  int NumberOfBins=sizeof(bins)/sizeof(bins[0])-1;
  TString InnerBarrelCut="fabs(ProbeEta)<0.8";
  //TString InnerBarrelCut="fabs(ProbeEta)>1 && fabs(ProbeEta)<1.5";
  TString OuterBarrelCut="fabs(ProbeEta)>0.8 && fabs(ProbeEta)<1.44";
  TString EndcapCut="fabs(ProbeEta)>1.57 && fabs(ProbeEta)<2.5";
  TFile *file_data = TFile::Open("TnPpairs_DATA.root");
  TTree *tree_data = (TTree*)file_data->Get("Events");
  tree_data->SetBranchStatus("*",0);
  tree_data->SetBranchStatus("Diele_mass",1);
  tree_data->SetBranchStatus("ele_pt",1);
  tree_data->SetBranchStatus("ProbePt",1);
  tree_data->SetBranchStatus("ProbeEta",1);
  tree_data->SetBranchStatus("ProbePass",1);
  gROOT->cd();

  TFile *file_MC = TFile::Open("TnPpairs_MC.root");
  TTree *tree_MC = (TTree*)file_MC->Get("Events");
  tree_MC->SetBranchStatus("*",0);
  tree_MC->SetBranchStatus("Diele_mass",1);
  tree_MC->SetBranchStatus("ele_pt",1);
  tree_MC->SetBranchStatus("ProbePt",1);
  tree_MC->SetBranchStatus("ProbeEta",1);
  tree_MC->SetBranchStatus("ProbePass",1);
  gROOT->cd();

  //Passing and failling probes (data and MC)
  TTree* PassingProbesData[NumberOfBins];
  TTree* FaillingProbesData[NumberOfBins];
  TTree* PassingProbesMC[NumberOfBins];
  TTree* FaillingProbesMC[NumberOfBins];

  if(FitInnerBarrel)
  for(i=0;i<NumberOfBins;i++)
  {
     PassingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==1 && " + InnerBarrelCut,bins[i],bins[i+1]));
     FaillingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==0 && " + InnerBarrelCut,bins[i],bins[i+1]));
     PassingProbesMC[i] = tree_MC->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==1 && " + InnerBarrelCut,bins[i],bins[i+1]));
     FaillingProbesMC[i] = tree_MC->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==0 && " + InnerBarrelCut,bins[i],bins[i+1]));
  }

  if(FitOuterBarrel)
  for(i=0;i<NumberOfBins;i++)
  {
     PassingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==1 && " + OuterBarrelCut,bins[i],bins[i+1]));
     FaillingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==0 && " + OuterBarrelCut,bins[i],bins[i+1]));
     PassingProbesMC[i] = tree_MC->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==1 && " + OuterBarrelCut,bins[i],bins[i+1]));
     FaillingProbesMC[i] = tree_MC->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==0 && " + OuterBarrelCut,bins[i],bins[i+1]));
  }

  if(FitEndcap)
  for(i=0;i<NumberOfBins;i++)
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

  //Efficiency(YieldErrorDataPass,YieldErrorDataFail,YieldErrorMCPass, YieldErrorMCFail, bins, NumberOfBins);

  //std::cout<<Efficiency(YieldErrorMCPass, YieldErrorMCFail, bins, NumberOfBins)->GetEfficiency(2)<<std::endl;

  TEfficiency* Eff_MC = Efficiency(YieldErrorMCPass, YieldErrorMCFail, bins, NumberOfBins);
  TEfficiency* Eff_DATA = Efficiency(YieldErrorDataPass, YieldErrorDataFail, bins, NumberOfBins);
  for(i=1;i<=3;i++)
  {
    std::cout<<Eff_MC->GetEfficiency(i)<<" "<<Eff_MC->GetEfficiencyErrorLow(i)<<" "<<Eff_MC->GetEfficiencyErrorUp(i)<<std::endl;
  }
  for(i=1;i<=3;i++)
  {
    std::cout<<Eff_DATA->GetEfficiency(i)<<" "<<Eff_DATA->GetEfficiencyErrorLow(i)<<" "<<Eff_DATA->GetEfficiencyErrorUp(i)<<std::endl;
  }
  std::cout<<YieldErrorMCPass.at(0).first<<" "<<YieldErrorMCPass.at(0).second<<std::endl;
  std::cout<<YieldErrorMCFail.at(0).first<<" "<<YieldErrorMCFail.at(0).second<<std::endl;
  auto legend = new TLegend(0.9,0.9,1.0,1.0);
  auto mg  = new TMultiGraph();
  Eff_MC->SetLineColor(kRed);
  Eff_DATA->SetLineColor(kBlue);
  TCanvas* oi = new TCanvas();
  oi->cd();
  Eff_MC->Draw();
  gPad->Update();
  auto graph_MC = Eff_MC->GetPaintedGraph();
  mg->Add(graph_MC);
  Eff_DATA->Draw();
  gPad->Update();
  auto graph_data = Eff_DATA->GetPaintedGraph();
  mg->Add(graph_data);
  mg->Draw("ap");
  mg->GetXaxis()->SetTitle("pT/GeV");
  mg->GetYaxis()->SetTitle("Efficiency");
  legend->AddEntry(graph_data,"data");
  legend->AddEntry(graph_MC,"MC");
  legend->Draw();
  oi->SaveAs("plots/Efficiency.pdf");
}
