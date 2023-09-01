#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "RooCrystalBall.cxx"
#include "RooCrystalBall.h"
#include "MakeFitData.cpp"
#include "MakeFitDataAltBgk.cpp"
#include "MakeFitDataAltSig.cpp"
#include "MakeFitDataNom.cpp"
#include "MakeFitMC.cpp"
#include "GetFitMC.cpp"
#include "GetFitData.cpp"
#include "Efficiency.cpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <stdlib.h>


void TnPAnalyzer()
{
  int i;
  //std::vector<std::pair<float,float>> YieldErrorDataPass, YieldErrorDataFail, YieldErrorMCPass, YieldErrorMCFail;

  std::vector<std::pair<float,float>> YieldErrorDataPass, YieldErrorDataFail, YieldErrorMCPass, YieldErrorMCFail;
  std::vector<std::pair<float,float>> YieldErrorDataPassBkg, YieldErrorDataFailBkg, YieldErrorDataPassSig, YieldErrorDataFailSig;
  std::vector<std::pair<float,float>> YieldErrorDataPassTag, YieldErrorDataFailTag, YieldErrorMCPassTag, YieldErrorMCFailTag;
  
  float bins[]={5,7,10,15};  //this part is configurable, select bins for fitting
  bool FitInnerBarrel=1, FitOuterBarrel=0, FitEndcap=0;
  int NumberOfBins=sizeof(bins)/sizeof(bins[0])-1;
  float mvaCutValue=0.97; //altTag
  int etaBin=0;  
  TString altCutData=TString::Format(" TagMVA>%f ",mvaCutValue);
  TString altCutMC=TString::Format(" TagMVA>%f && TagPt>7 ",mvaCutValue);
  TString etaBinName = "";
  TString EtaCut="";

  std::ofstream Eff_file("plots/Eff_output.txt"); // defining .txt output
  std::ofstream SF_file("plots/SF_output.txt"); // defining .txt output

  Eff_file<<"etaBin\tptBin\td_eff_n\td_eff_n_err\tm_eff\tm_eff_err\td_eff_b\td_eff_b_err\td_eff_s\td_eff_s_err\td_eff_t\td_eff_t_err\tm_eff_t\tm_eff_t_err"<<std::endl;
  SF_file<<"etaBin\tptBin\tSF_n\tSF_n_err\tSF_b\tSF_b_err\tSF_s\tSF_s_err\tSF_t\tSF_t_err\tSF_n\tSF_err_tot"<<std::endl;
  

  TFile *file_data = TFile::Open("/eos/user/a/angaile/bParkingRootFiles/wpHZZSum18/TnPpairs_UL_Data_wpHZZSum.root");
  //TFile *file_data = TFile::Open("TnPpairs_DATA.root");
  TTree *tree_data = (TTree*)file_data->Get("Events");
  tree_data->SetBranchStatus("*",0);
  tree_data->SetBranchStatus("Diele_mass",1);
  //tree_data->SetBranchStatus("ele_ip3D_match",1);
  tree_data->SetBranchStatus("ele_pt",1);
  tree_data->SetBranchStatus("ProbePt",1);
  //tree_data->SetBranchStatus("ProbeDxy",1);
  tree_data->SetBranchStatus("ProbeEta",1);
  tree_data->SetBranchStatus("ProbePass",1);
  tree_data->SetBranchStatus("ProbeRelISO",1);
  tree_data->SetBranchStatus("TagRelISO",1);
  tree_data->SetBranchStatus("TagMVA",1);
  tree_data->SetBranchStatus("TagPt",1);
  gROOT->cd();

  TFile *file_MC = TFile::Open("/eos/user/a/angaile/bParkingRootFiles/wpHZZSum18/TnPpairs_UL_MC_wpHZZSum.root");
  //TFile *file_MC = TFile::Open("TnPpairs_MC.root");
  TTree *tree_MC = (TTree*)file_MC->Get("Events");
  tree_MC->SetBranchStatus("*",0);
  tree_MC->SetBranchStatus("Diele_mass",1);
  tree_MC->SetBranchStatus("ele_pt",1);
  //tree_MC->SetBranchStatus("ele_ip3D_match",1);
  tree_MC->SetBranchStatus("ProbePt",1);
  //tree_MC->SetBranchStatus("ProbeDxy",1);
  tree_MC->SetBranchStatus("ProbeEta",1);
  tree_MC->SetBranchStatus("ProbePass",1);
  tree_MC->SetBranchStatus("ProbeRelISO",1);
  tree_MC->SetBranchStatus("TagRelISO",1);
  tree_MC->SetBranchStatus("TagMVA",1);
  tree_MC->SetBranchStatus("TagPt",1);
  gROOT->cd();

  for(etaBin=1;etaBin<=3;etaBin++)
	{// rolling through etas
    if(etaBin==1){
      std::cout << "\nPlotting Inner barrel\n "  << std::endl;
      FitInnerBarrel=1;
      FitOuterBarrel=0;
      FitEndcap=0;
      etaBinName="inner";
      EtaCut=" fabs(ProbeEta)<0.8 ";
    } else if(etaBin==2){
      std::cout << "\nPlotting Outter barrel\n "  << std::endl;
      FitInnerBarrel=0;
      FitOuterBarrel=1;
      FitEndcap=0;
      etaBinName="outter";
      EtaCut=" fabs(ProbeEta)>0.8 && fabs(ProbeEta)<1.44 ";
    } else {
      std::cout << "\nPlotting Endcaps\n "  << std::endl;
      FitInnerBarrel=0;
      FitOuterBarrel=0;
      FitEndcap=1;
      etaBinName="endcaps";
      EtaCut=" fabs(ProbeEta)>1.57 && fabs(ProbeEta)<2.5 ";
    }

    //Passing and failling probes (data and MC)
    TTree* PassingProbesData[NumberOfBins];
    TTree* FaillingProbesData[NumberOfBins];
    TTree* PassingProbesMC[NumberOfBins];
    TTree* FaillingProbesMC[NumberOfBins];

    // Passing and failling probes for alt Tag (data and MC)
    TTree* tPassingProbesData[NumberOfBins];
    TTree* tFaillingProbesData[NumberOfBins];
    TTree* tPassingProbesMC[NumberOfBins];
    TTree* tFaillingProbesMC[NumberOfBins];

    string namesBins[NumberOfBins];

    for(i=0;i<NumberOfBins;i++)
    {
      PassingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==1 && " + EtaCut,bins[i],bins[i+1]));
      FaillingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==0 && " + EtaCut,bins[i],bins[i+1]));
      PassingProbesMC[i] = tree_MC->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==1 && " + EtaCut,bins[i],bins[i+1]));
      FaillingProbesMC[i] = tree_MC->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==0 && " + EtaCut,bins[i],bins[i+1]));

      tPassingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==1 && " + EtaCut +" && " + altCutData,bins[i],bins[i+1]));
      tFaillingProbesData[i] = tree_data->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==0 && " + EtaCut +" && " + altCutData,bins[i],bins[i+1]));
      tPassingProbesMC[i] = tree_MC->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==1 && " + EtaCut +" && " + altCutMC,bins[i],bins[i+1]));
      tFaillingProbesMC[i] = tree_MC->CopyTree(TString::Format("ProbePt>%f && ProbePt<%f && ProbePass==0 && " + EtaCut +" && " + altCutMC,bins[i],bins[i+1]));
    }


    //Fit for passing and failling probes in each bin
    for(i=0;i<NumberOfBins;i++)
    {
      std::pair<float,float> InfoMCPass, InfoMCFail, InfoMCPassT, InfoMCFailT;
      std::tuple<float,float,float,float,float,float,float,float,float> ParDataPass, ParDataFail,  ParDataPassT, ParDataFailT;

      InfoMCPass=GetFitMC(PassingProbesMC[i]);
      ParDataPass=GetFitData(PassingProbesData[i],InfoMCPass);
      InfoMCFail=GetFitMC(FaillingProbesMC[i]);
      ParDataFail=GetFitData(FaillingProbesData[i],InfoMCFail);

      InfoMCPassT=GetFitMC(tPassingProbesMC[i]);
      ParDataPassT=GetFitData(tPassingProbesData[i],InfoMCPassT);
      InfoMCFailT=GetFitMC(tFaillingProbesMC[i]);
      ParDataFailT=GetFitData(tFaillingProbesData[i],InfoMCFailT);


      // YieldErrorDataPass.push_back(MakeFitData(PassingProbesData[i], 1,  "plots/PassingProbesData_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf"));
      // YieldErrorDataPass.push_back(MakeFitData(PassingProbesData[i], 1,  "plots/PassingProbesData_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".root"));
      // YieldErrorDataFail.push_back(MakeFitData(FaillingProbesData[i], 1, "plots/FaillingProbesData_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf"));
      // YieldErrorDataFail.push_back(MakeFitData(FaillingProbesData[i], 1, "plots/FaillingProbesData_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".root"));
      // YieldErrorMCPass.push_back(MakeFitMC(PassingProbesMC[i], 1, "plots/PassingProbesMC_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf"));
      // YieldErrorMCPass.push_back(MakeFitMC(PassingProbesMC[i], 1, "plots/PassingProbesMC_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".root"));
      // YieldErrorMCFail.push_back(MakeFitMC(FaillingProbesMC[i], 1, "plots/FaillingProbesMC_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf"));
      // YieldErrorMCFail.push_back(MakeFitMC(FaillingProbesMC[i], 1, "plots/FaillingProbesMC_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".root"));

      // regular fits
      YieldErrorDataPass.push_back(MakeFitDataNom(PassingProbesData[i],ParDataPass, true, "plots/PassingProbesData_nom_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf",false));
      //YieldErrorDataPass.push_back(MakeFitDataNom(PassingProbesData[i],ParDataPass, true, "plots/PassingProbesData_nom_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".root",false));
      YieldErrorDataFail.push_back(MakeFitDataNom(FaillingProbesData[i],ParDataFail, true, "plots/FaillingProbesData_nom_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf",false));
      //YieldErrorDataFail.push_back(MakeFitDataNom(FaillingProbesData[i],ParDataFail, true, "plots/FaillingProbesData_nom_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".root",false));
      
      YieldErrorMCPass.push_back(MakeFitMC(PassingProbesMC[i], true, "plots/PassingProbesMC_nom_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf"));
      //YieldErrorMCPass.push_back(MakeFitMC(PassingProbesMC[i], true, "plots/PassingProbesMC_nom_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".root"));
      YieldErrorMCFail.push_back(MakeFitMC(FaillingProbesMC[i], true, "plots/FaillingProbesMC_nom_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf"));
      //YieldErrorMCFail.push_back(MakeFitMC(FaillingProbesMC[i], true, "plots/FaillingProbesMC_nom_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".root"));
    
      // alt tag fits
      YieldErrorDataPassTag.push_back(MakeFitDataNom(tPassingProbesData[i],ParDataPassT, true,  "plots/PassingProbesData_altT_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf",false));
      //YieldErrorDataPassTag.push_back(MakeFitDataNom(tPassingProbesData[i],ParDataPassT, true,  "plots/PassingProbesData_altT_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".root",false));
      YieldErrorDataFailTag.push_back(MakeFitDataNom(tFaillingProbesData[i],ParDataFailT, true, "plots/FaillingProbesData_altT_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf",false));
      //YieldErrorDataFailTag.push_back(MakeFitDataNom(tFaillingProbesData[i],ParDataFailT, true, "plots/FaillingProbesData_altT_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".root",false));
      YieldErrorMCPassTag.push_back(MakeFitMC(tPassingProbesMC[i], true, "plots/PassingProbesMC_altT_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf"));
      //YieldErrorMCPassTag.push_back(MakeFitMC(tPassingProbesMC[i], true, "plots/PassingProbesMC_altT_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".root"));
      YieldErrorMCFailTag.push_back(MakeFitMC(tFaillingProbesMC[i], true, "plots/FaillingProbesMC_altT_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf"));
      //YieldErrorMCFailTag.push_back(MakeFitMC(tFaillingProbesMC[i], true, "plots/FaillingProbesMC_altT_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".root"));
    
      // alt bgk fits
      YieldErrorDataPassBkg.push_back(MakeFitDataAltBgk(PassingProbesData[i],ParDataPass, true,  "plots/PassingProbesData_altB_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf",false));
      //YieldErrorDataPassBkg.push_back(MakeFitDataAltBgk(PassingProbesData[i],ParDataPass, true,  "plots/PassingProbesData_altB_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".root",false));
      YieldErrorDataFailBkg.push_back(MakeFitDataAltBgk(FaillingProbesData[i],ParDataFail, true, "plots/FaillingProbesData_altB_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf",false));
      //YieldErrorDataFailBkg.push_back(MakeFitDataAltBgk(FaillingProbesData[i],ParDataFail, true, "plots/FaillingProbesData_altB_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".root",false));

      // alt sign fits
      //std::cout << "\n\n\nStarting alt\n\n\n "  << std::endl;
      YieldErrorDataPassSig.push_back(MakeFitDataAltSig(PassingProbesData[i],ParDataPass, true,  "plots/PassingProbesData_altS_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf",false));
      //YieldErrorDataPassSig.push_back(MakeFitDataAltSig(PassingProbesData[i],ParDataPass, true,  "plots/PassingProbesData_altS_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".root",false));
      YieldErrorDataFailSig.push_back(MakeFitDataAltSig(FaillingProbesData[i],ParDataFail, true, "plots/FaillingProbesData_altS_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".pdf",false));
      //YieldErrorDataFailSig.push_back(MakeFitDataAltSig(FaillingProbesData[i],ParDataFail, true, "plots/FaillingProbesData_altS_"+std::to_string(etaBin)+"_"+std::to_string((int)bins[i])+"_"+std::to_string((int)bins[i+1])+".root",false));


    }

    //TEfficiency* Eff_MC = Efficiency(YieldErrorMCPass, YieldErrorMCFail, bins, NumberOfBins);
    //TEfficiency* Eff_DATA = Efficiency(YieldErrorDataPass, YieldErrorDataFail, bins, NumberOfBins);

    TEfficiency* Eff_MC = Efficiency(YieldErrorMCPass, YieldErrorMCFail, bins, NumberOfBins);
    TEfficiency* Eff_MC_altTag = Efficiency(YieldErrorMCPassTag, YieldErrorMCFailTag, bins, NumberOfBins);

    TEfficiency* Eff_DATA = Efficiency(YieldErrorDataPass, YieldErrorDataFail, bins, NumberOfBins);
    TEfficiency* Eff_DATA_altBkg = Efficiency(YieldErrorDataPassBkg, YieldErrorDataFailBkg, bins, NumberOfBins);
    TEfficiency* Eff_DATA_altSig = Efficiency(YieldErrorDataPassSig, YieldErrorDataFailSig, bins, NumberOfBins);
    TEfficiency* Eff_DATA_altTag = Efficiency(YieldErrorDataPassTag, YieldErrorDataFailTag, bins, NumberOfBins);
    
    
    for(i=1;i<=NumberOfBins;i++) // i<=3 (i=1;i<=3;i++)
    {
      std::cout<<Eff_MC->GetEfficiency(i)<<" "<<Eff_MC->GetEfficiencyErrorLow(i)<<" "<<Eff_MC->GetEfficiencyErrorUp(i)<<std::endl;
      
      Eff_file << etaBin << "\t" << i << "\t" << Eff_DATA->GetEfficiency(i) << "\t" <<Eff_DATA->GetEfficiencyErrorUp(i)<< "\t" <<Eff_MC->GetEfficiency(i)<< "\t" <<Eff_MC->GetEfficiencyErrorUp(i)<< "\t" << Eff_DATA_altBkg->GetEfficiency(i) << "\t" <<Eff_DATA_altBkg->GetEfficiencyErrorUp(i)<< "\t" << Eff_DATA_altSig->GetEfficiency(i) << "\t" <<Eff_DATA_altSig->GetEfficiencyErrorUp(i)<< "\t" << Eff_DATA_altTag->GetEfficiency(i) << "\t" <<Eff_DATA_altTag->GetEfficiencyErrorUp(i)<< "\t" << Eff_MC_altTag->GetEfficiency(i) << "\t" <<Eff_MC_altTag->GetEfficiencyErrorUp(i)<< std::endl;
    }

    //std::cout<<YieldErrorMCPass.at(0).first<<" "<<YieldErrorMCPass.at(0).second<<std::endl;
    //std::cout<<YieldErrorMCFail.at(0).first<<" "<<YieldErrorMCFail.at(0).second<<std::endl;

    //scale factors
    TH1F* SF = new TH1F("scale factor", "scale factor",NumberOfBins, bins);
    for(i=1;i<=NumberOfBins;i++)
    {
      
      double SF_n = Eff_DATA->GetEfficiency(i)*1.0/Eff_MC->GetEfficiency(i);
      double SF_altB = Eff_DATA_altBkg->GetEfficiency(i)*1.0/Eff_MC->GetEfficiency(i);
      double SF_altS = Eff_DATA_altSig->GetEfficiency(i)*1.0/Eff_MC->GetEfficiency(i);
      double SF_altT = Eff_DATA_altTag->GetEfficiency(i)*1.0/Eff_MC_altTag->GetEfficiency(i);

      double SF_n_err = (Eff_DATA->GetEfficiency(i)*1.0/Eff_MC->GetEfficiency(i))*sqrt(pow((Eff_DATA->GetEfficiencyErrorUp(i)/Eff_DATA->GetEfficiency(i)),2)+pow((Eff_MC->GetEfficiencyErrorUp(i)/Eff_MC->GetEfficiency(i)),2));
      double SF_altB_err = (Eff_DATA_altBkg->GetEfficiency(i)*1.0/Eff_MC->GetEfficiency(i))*sqrt(pow((Eff_DATA_altBkg->GetEfficiencyErrorUp(i)/Eff_DATA_altBkg->GetEfficiency(i)),2)+pow((Eff_MC->GetEfficiencyErrorUp(i)/Eff_MC->GetEfficiency(i)),2));
      double SF_altS_err = (Eff_DATA_altSig->GetEfficiency(i)*1.0/Eff_MC->GetEfficiency(i))*sqrt(pow((Eff_DATA_altSig->GetEfficiencyErrorUp(i)/Eff_DATA_altSig->GetEfficiency(i)),2)+pow((Eff_MC->GetEfficiencyErrorUp(i)/Eff_MC->GetEfficiency(i)),2));
      double SF_altT_err = (Eff_DATA_altTag->GetEfficiency(i)*1.0/Eff_MC_altTag->GetEfficiency(i))*sqrt(pow((Eff_DATA_altTag->GetEfficiencyErrorUp(i)/Eff_DATA_altTag->GetEfficiency(i)),2)+pow((Eff_MC_altTag->GetEfficiencyErrorUp(i)/Eff_MC_altTag->GetEfficiency(i)),2));
      
      double SF_n_err_total = sqrt(pow(SF_n_err,2)+pow(SF_n-SF_altB,2)+pow(SF_n-SF_altS,2)+pow(SF_n-SF_altT,2));

      SF->SetBinContent(i,Eff_DATA->GetEfficiency(i)*1.0/Eff_MC->GetEfficiency(i));
      //double SF_error=sqrt(pow(Eff_MC->GetEfficiency(i)*Eff_DATA->GetEfficiencyErrorUp(i),2)+pow(Eff_DATA->GetEfficiency(i)*Eff_MC->GetEfficiencyErrorUp(i),2))/pow(Eff_MC->GetEfficiency(i),2);
      
      SF->SetBinError(i, SF_n_err_total);
      // My extra line for SF
      //std::cout << "Bin value " << i << ": " << Eff_DATA->GetEfficiency(i)*1.0/Eff_MC->GetEfficiency(i) << " | Error : " << SF_error << std::endl;
      SF_file << etaBin << "\t" << i << "\t" << SF_n << "\t" << SF_n_err << "\t" << SF_altB << "\t" << SF_altB_err << "\t" << SF_altS << "\t" << SF_altS_err << "\t" << SF_altT << "\t" << SF_altT_err << "\t" << SF_n << "\t" << SF_n_err_total << std::endl;
      //
    }
    
    Eff_file.close();
    SF_file.close();

    SF->GetXaxis()->SetTitle("pT/GeV");
    SF->GetYaxis()->SetTitle("SF");
    TCanvas* SF_canvas = new TCanvas();
    gStyle->SetOptStat(0);
    SF->Draw("E");
    SF_canvas->SaveAs("plots/ScaleFactor.pdf");

    //efficiency plot
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

    //Save eff object to root file
    TFile* eff = TFile::Open("efficiency.root","RECREATE");
    eff->WriteObject(Eff_MC, "Eff_MC");
    eff->WriteObject(Eff_DATA, "Eff_DATA");

    //Read Tefficiency object from root file
    TEfficiency* tmp = (TEfficiency*)eff->Get("Eff_MC");
    std::cout<<tmp->GetEfficiency(1)<<std::endl;
  } 
}