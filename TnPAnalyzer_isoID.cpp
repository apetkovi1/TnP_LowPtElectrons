#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "RooCrystalBall.cxx"
#include "RooCrystalBall.h"
#include "MakeFitData.cpp"
#include "MakeFitDataAltBgk.cpp"
#include "MakeFitDataAltSig.cpp"
#include "MakeFitMC.cpp"
#include "GetFitMC.cpp"
#include "GetFitData.cpp"
#include "MakeFitDataNom.cpp"
#include "Efficiency.cpp"
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <tuple>
#include <vector>

void TnPAnalyzerGrandUL_wpHZZSum()
{
  int i;
  string gribuTag="tagNoIso"; //decoration for png names
  string folderName="HZZtest";//idChoice;
  bool saglabat=1;
  bool doEfficiency=1;
  string outFileName="plots/"+folderName+"/UL_altRez_eff.txt";
  // creating output file to store the results
  std::ofstream outfile;
  outfile.open(outFileName, std::ios_base::app); // append instead of overwrite
  
  outfile<<"eta\tpt\trelIso\td_eff_n\td_eff_n_err\tm_eff\tm_eff_err\td_eff_b\td_eff_b_err\td_eff_s\td_eff_s_err\td_eff_t\td_eff_t_err\tm_eff_t\tm_eff_t_err\tcomment"<<std::endl;
  float MassConstrainPass, MassConstrainFail;  
  float mvaCutValue=0.97;//altTag  
  TString altCutData=TString::Format(" TagMVA>%f ",mvaCutValue);
  TString altCutMC=TString::Format(" TagMVA>%f && TagPt>7 ",mvaCutValue);
  
  float ptBins[]={5,7,10,15};  //this part is configurable, select bins for fitting
  bool FitInnerBarrel=1, FitOuterBarrel=0, FitEndcap=0;
  TString EtaCut="";
  TString isoCut="";
  TString ptCut="";
  int etaIzvele=0;
  int isoIzvele=0;
  int ptIzvele=0;
  TString gribu="";
  
  
  // importing DATA file
  TFile *file_data = TFile::Open("/eos/user/a/angaile/bParkingRootFiles/wpHZZSum18/TnPpairs_UL_Data_wpHZZSum.root");
  TTree *tree_data = (TTree*)file_data->Get("Events");
  tree_data->SetBranchStatus("*",0);
  tree_data->SetBranchStatus("Diele_mass",1);
  tree_data->SetBranchStatus("ele_pt",1);
  tree_data->SetBranchStatus("ProbePt",1);
  tree_data->SetBranchStatus("ProbeEta",1);
  tree_data->SetBranchStatus("ProbePass",1);
  tree_data->SetBranchStatus("ProbeRelISO",1);
  tree_data->SetBranchStatus("TagMVA",1);
  gROOT->cd();

  // importing MC file
  TFile *file_MC = TFile::Open("/eos/user/a/angaile/bParkingRootFiles/wpHZZSum18/TnPpairs_UL_MC_wpHZZSum.root");
  //TFile *file_MC = TFile::Open("TnPpairs_UL_MC_pIso_tNoIso.root");
  TTree *tree_MC = (TTree*)file_MC->Get("Events");
  tree_MC->SetBranchStatus("*",0);
  tree_MC->SetBranchStatus("Diele_mass",1);
  tree_MC->SetBranchStatus("ele_pt",1);
  tree_MC->SetBranchStatus("ProbePt",1);
  tree_MC->SetBranchStatus("TagPt",1);
  tree_MC->SetBranchStatus("ProbeEta",1);
  tree_MC->SetBranchStatus("ProbePass",1);
  tree_MC->SetBranchStatus("ProbeRelISO",1);
  //  tree_MC->SetBranchStatus("TagMVAIso",1);
  tree_MC->SetBranchStatus("TagMVA",1);
  gROOT->cd();
  
  for(etaIzvele=1;etaIzvele<=3;etaIzvele++)
	{// rolling through etas
    if(etaIzvele==1){
      std::cout << "\nPlotting Inner barrel\n "  << std::endl;
      FitInnerBarrel=1;
      FitOuterBarrel=0;
      FitEndcap=0;
      gribu="inner";
      EtaCut=" fabs(ProbeEta)<0.8 ";
    } else if(etaIzvele==2){
      std::cout << "\nPlotting Outter barrel\n "  << std::endl;
      FitInnerBarrel=0;
      FitOuterBarrel=1;
      FitEndcap=0;
      gribu="outter";
      EtaCut=" fabs(ProbeEta)>0.8 && fabs(ProbeEta)<1.44 ";
    } else {
      std::cout << "\nPlotting Endcaps\n "  << std::endl;
      FitInnerBarrel=0;
      FitOuterBarrel=0;
      FitEndcap=1;
      gribu="endcaps";
      EtaCut=" fabs(ProbeEta)>1.57 && fabs(ProbeEta)<2.5 ";
    }
    
    std::cout <<EtaCut << std::endl; 
    int NumberOfPtBins=sizeof(ptBins)/sizeof(ptBins[0])-1;
    
    for(ptIzvele=0;ptIzvele<NumberOfPtBins;ptIzvele++)
    {  
      ptCut=TString::Format(" ProbePt>%f && ProbePt<%f ",ptBins[ptIzvele],ptBins[ptIzvele+1]);
      // for each eta region
      std::vector<std::pair<float,float>> YieldErrorDataPass, YieldErrorDataFail, YieldErrorMCPass, YieldErrorMCFail;
      std::vector<std::pair<float,float>> YieldErrorDataPassBkg, YieldErrorDataFailBkg, YieldErrorDataPassSig, YieldErrorDataFailSig;
      std::vector<std::pair<float,float>> YieldErrorDataPassTag, YieldErrorDataFailTag, YieldErrorMCPassTag, YieldErrorMCFailTag;


      float isoBins[5]={0};
      if(ptBins[ptIzvele+1]== 7 )
      {
        if(etaIzvele==1) {
          float pisoBins[]={0,0.001,0.155,0.330,0.5};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
        if(etaIzvele==2) {
          float pisoBins[]={0,0.001,0.155,0.324,0.5};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
        if(etaIzvele==3) {
          float pisoBins[]={0,0.001,0.204,0.5,2};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
      } else if(ptBins[ptIzvele+1]== 10 )
      {
        if(etaIzvele==1) {
          float pisoBins[]={0,0.081,0.277,0.452,0.6};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
        if(etaIzvele==2) {
          float pisoBins[]={0,0.095,0.286,0.455,0.6};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
        if(etaIzvele==3) {
          float pisoBins[]={0,0.065,0.293,0.6,2};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
      } else if(ptBins[ptIzvele+1]== 15 )
      {
        if(etaIzvele==1) {
          float pisoBins[]={0,0.040,0.141,0.226,0.3};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
        if(etaIzvele==2) {
          float pisoBins[]={0,0.037,0.139,0.224,0.3};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
        if(etaIzvele==3) {
          float pisoBins[]={0,0.040,0.164,0.3,2};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
      }

      int NumberOfIsoBins=sizeof(isoBins)/sizeof(isoBins[0])-1;
      int NumberOfBins=NumberOfIsoBins;//NumberOfPtBins*NumberOfIsoBins;

      // Passing and failling probes (data and MC)
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
      string namesPt="pt_"+std::to_string((int)ptBins[ptIzvele])+"_"+std::to_string((int)ptBins[ptIzvele+1]);
      i=0;
      for(isoIzvele=0;isoIzvele<NumberOfIsoBins;isoIzvele++)
      {// rolling through iso
        isoCut=TString::Format(" ProbeRelISO>=%f && ProbeRelISO<%f ",isoBins[isoIzvele],isoBins[isoIzvele+1]);
      
        namesBins[i] = namesPt+"_relIso_"+std::to_string((int)(isoBins[isoIzvele]*100))+"_"+std::to_string((int)(isoBins[isoIzvele+1]*100));
      
        PassingProbesData[i] = tree_data->CopyTree(ptCut+"&&"+isoCut + " && ProbePass==1 && " + EtaCut);
        FaillingProbesData[i] = tree_data->CopyTree(ptCut+"&&"+isoCut + " && ProbePass==0 && " + EtaCut);
        
        PassingProbesMC[i] = tree_MC->CopyTree(ptCut+"&&"+isoCut + " && ProbePass==1 && " + EtaCut + " && TagPt>7");
        FaillingProbesMC[i] = tree_MC->CopyTree(ptCut+"&&"+isoCut + " && ProbePass==0 && " + EtaCut + " && TagPt>7");
         
        tPassingProbesData[i] = tree_data->CopyTree(ptCut+"&&"+isoCut + " && ProbePass==1 && " + EtaCut+" && "+ altCutData); 
        tFaillingProbesData[i] = tree_data->CopyTree(ptCut+"&&"+isoCut + " && ProbePass==0 && " + EtaCut+" && "+ altCutData);
        tPassingProbesMC[i] = tree_MC->CopyTree(ptCut+"&&"+isoCut + " && ProbePass==1 && " + EtaCut+" && "+ altCutMC);
        tFaillingProbesMC[i] = tree_MC->CopyTree(ptCut+"&&"+isoCut + " && ProbePass==0 && "+ EtaCut+" && "+ altCutMC);
        i++;
      }

      //Fit for passing and failling probes in each bin
      //std::vector<std::pair<float,float>> InfoMCPass, InfoMCFail;
      //std::tuple<float,float,float,float,float,float,float,float,float> ParDataPass, ParDataFail;
      //for(i=0;i<1;i++)
      for(i=0;i<NumberOfBins;i++)//(i=0;i<NumberOfBins;i++)
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

      
        string plotNamePart="plots/"+folderName+"/UL_eta_"+std::to_string(etaIzvele)+"_"+namesBins[i]+"_"+gribuTag;
        if(isoBins[i+1]==2)
          saglabat=0;
        else
          saglabat=1;
        // regular fits
        YieldErrorDataPass.push_back(MakeFitDataNom(PassingProbesData[i],ParDataPass, saglabat, plotNamePart+"_nom_PassingProbesData.png",false));
        //std::cout << "it was passing data " << namesBins[i] << "\n\n\n\n\n";
        YieldErrorDataFail.push_back(MakeFitDataNom(FaillingProbesData[i],ParDataFail, saglabat, plotNamePart+"_nom_FaillingProbesData.png",false));
        //std::cout << "it was failing data " << namesBins[i] << "\n\n\n\n\n";    
        YieldErrorMCPass.push_back(MakeFitMC(PassingProbesMC[i], saglabat, plotNamePart+"_nom_PassingProbesMC.png"));
        YieldErrorMCFail.push_back(MakeFitMC(FaillingProbesMC[i], saglabat, plotNamePart+"_nom_FaillingProbesMC.png"));
      
        // alt tag fits
        YieldErrorDataPassTag.push_back(MakeFitDataNom(tPassingProbesData[i],ParDataPassT, saglabat,  plotNamePart+"_altTag_PassingProbesData.png",false));
        YieldErrorDataFailTag.push_back(MakeFitDataNom(tFaillingProbesData[i],ParDataFailT, saglabat, plotNamePart+"_altTag_FaillingProbesData.png",false));
        YieldErrorMCPassTag.push_back(MakeFitMC(tPassingProbesMC[i], saglabat, plotNamePart+"_altTag_PassingProbesMC.png"));
        YieldErrorMCFailTag.push_back(MakeFitMC(tFaillingProbesMC[i], saglabat, plotNamePart+"_altTag_FaillingProbesMC.png"));
      
        // alt bgk fits
        YieldErrorDataPassBkg.push_back(MakeFitDataAltBgk(PassingProbesData[i],ParDataPass, saglabat,  plotNamePart+"_altBgk_PassingProbesData.png",false));
        YieldErrorDataFailBkg.push_back(MakeFitDataAltBgk(FaillingProbesData[i],ParDataFail, saglabat, plotNamePart+"_altBgk_FaillingProbesData.png",false));
        
        // alt sign fits
        //std::cout << "\n\n\nStarting alt\n\n\n "  << std::endl;
        YieldErrorDataPassSig.push_back(MakeFitDataAltSig(PassingProbesData[i],ParDataPass, saglabat,  plotNamePart+"_altSig_PassingProbesData.png",false));
        YieldErrorDataFailSig.push_back(MakeFitDataAltSig(FaillingProbesData[i],ParDataFail, saglabat, plotNamePart+"_altSig_FaillingProbesData.png",false));
      }

      // here it should just be relIso
      TEfficiency* Eff_MC = Efficiency(YieldErrorMCPass, YieldErrorMCFail, isoBins, NumberOfBins);
      TEfficiency* Eff_MC_altTag = Efficiency(YieldErrorMCPassTag, YieldErrorMCFailTag, isoBins, NumberOfBins);

      TEfficiency* Eff_DATA = Efficiency(YieldErrorDataPass, YieldErrorDataFail, isoBins, NumberOfBins);
      // std::cout << "it was Eff calc " << "\n";  
      TEfficiency* Eff_DATA_altBkg = Efficiency(YieldErrorDataPassBkg, YieldErrorDataFailBkg, isoBins, NumberOfBins);
      TEfficiency* Eff_DATA_altSig = Efficiency(YieldErrorDataPassSig, YieldErrorDataFailSig, isoBins, NumberOfBins);
      TEfficiency* Eff_DATA_altTag = Efficiency(YieldErrorDataPassTag, YieldErrorDataFailTag, isoBins, NumberOfBins);


      
      for(i=1;i<=NumberOfIsoBins;i++)
      {
        if(isoBins[i+1]==2)
          continue;
        outfile<<std::to_string(etaIzvele)<<"\t"<<std::to_string((int)ptBins[ptIzvele+1])<<"\t"<<std::to_string((float)isoBins[i-1])<<"\t"<<Eff_DATA->GetEfficiency(i)<<"\t"<<Eff_DATA->GetEfficiencyErrorUp(i)<<"\t"<<Eff_MC->GetEfficiency(i)<<"\t"<<Eff_MC->GetEfficiencyErrorUp(i)<<"\t"<<Eff_DATA_altBkg->GetEfficiency(i)<<"\t"<<Eff_DATA_altBkg->GetEfficiencyErrorUp(i)<<"\t"<<Eff_DATA_altSig->GetEfficiency(i)<<"\t"<<Eff_DATA_altSig->GetEfficiencyErrorUp(i)<<"\t"<<Eff_DATA_altTag->GetEfficiency(i)<<"\t"<<Eff_DATA_altTag->GetEfficiencyErrorUp(i)<<"\t"<<Eff_MC_altTag->GetEfficiency(i)<<"\t"<<Eff_MC_altTag->GetEfficiencyErrorUp(i)<<"\t"<<namesBins[i-1]<<std::endl;
      }
      
      //Save eff object to root file
      //  TString effName="plots/"+idChoice+"/UL_eta_"+gribu+"_"+namesPt+"_efficiency_"+gribuTag+".root";
      TFile* eff = TFile::Open("plots/"+folderName+"/UL_eta_"+gribu+"_"+namesPt+"_efficiency_"+gribuTag+".root","RECREATE");
      eff->WriteObject(Eff_MC, "Eff_MC");
      eff->WriteObject(Eff_MC_altTag, "Eff_MC_altTag");
      eff->WriteObject(Eff_DATA, "Eff_DATA");
      eff->WriteObject(Eff_DATA_altSig, "Eff_DATA_altSig");
      eff->WriteObject(Eff_DATA_altBkg, "Eff_DATA_altBkg");
      eff->WriteObject(Eff_DATA_altTag, "Eff_DATA_altTag");
      

      if(doEfficiency){
        //efficiency plot regular fit
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
        mg->GetXaxis()->SetTitle("relIso");
        mg->GetYaxis()->SetTitle("Efficiency "+namesPt+", "+gribu);
        legend->AddEntry(graph_data,"data");
        legend->AddEntry(graph_MC,"MC");
        legend->Draw();
        oi->SaveAs("plots/"+folderName+"/UL_eta_"+gribu+"_"+namesPt+"_Efficiency_"+gribuTag+"_nom.png");

        //efficiency plot alternative background
        auto legendb = new TLegend(0.9,0.9,1.0,1.0);
        auto mgb  = new TMultiGraph();
        Eff_MC->SetLineColor(kRed);
        Eff_DATA_altBkg->SetLineColor(kBlue);
        TCanvas* oib = new TCanvas();
        oib->cd();
        Eff_MC->Draw();
        gPad->Update();
        auto graph_MCb = Eff_MC->GetPaintedGraph();
        mgb->Add(graph_MCb);
        Eff_DATA_altBkg->Draw();
        gPad->Update();
        auto graph_datab = Eff_DATA_altBkg->GetPaintedGraph();
        mgb->Add(graph_datab);
        mgb->Draw("ap");
        mgb->GetXaxis()->SetTitle("relIso");
        mgb->GetYaxis()->SetTitle("Efficiency_altBkg "+namesPt+", "+gribu);
        legendb->AddEntry(graph_datab,"data");
        legendb->AddEntry(graph_MCb,"MC");
        legendb->Draw();
        oib->SaveAs("plots/"+folderName+"/UL_eta_"+gribu+"_"+namesPt+"_Efficiency_"+gribuTag+"_altBkg.png");
        
        //efficiency plot alternative signal
        auto legends = new TLegend(0.9,0.9,1.0,1.0);
        auto mgs  = new TMultiGraph();
        Eff_MC->SetLineColor(kRed);
        Eff_DATA_altSig->SetLineColor(kBlue);
        TCanvas* ois = new TCanvas();
        ois->cd();
        Eff_MC->Draw();
        gPad->Update();
        auto graph_MCs = Eff_MC->GetPaintedGraph();
        mgs->Add(graph_MCs);
        Eff_DATA_altSig->Draw();
        gPad->Update();
        auto graph_datas = Eff_DATA_altSig->GetPaintedGraph();
        mgs->Add(graph_datas);
        mgs->Draw("ap");
        mgs->GetXaxis()->SetTitle("relIso");
        mgs->GetYaxis()->SetTitle("Efficiency_altSig "+namesPt+", "+gribu);
        legends->AddEntry(graph_datas,"data");
        legends->AddEntry(graph_MCs,"MC");
        legends->Draw();
        ois->SaveAs("plots/"+folderName+"/UL_eta_"+gribu+"_"+namesPt+"_Efficiency_"+gribuTag+"_altSig.png");
        
        //efficiency plot alternative tag
        auto legendt = new TLegend(0.9,0.9,1.0,1.0);
        auto mgt  = new TMultiGraph();
        Eff_MC_altTag->SetLineColor(kRed);
        Eff_DATA_altTag->SetLineColor(kBlue);
        TCanvas* oit = new TCanvas();
        oit->cd();
        Eff_MC_altTag->Draw();
        gPad->Update();
        auto graph_MC_altTag = Eff_MC_altTag->GetPaintedGraph();
        mgt->Add(graph_MC_altTag);
        Eff_DATA_altTag->Draw();
        gPad->Update();
        auto graph_data_altTag = Eff_DATA_altTag->GetPaintedGraph();
        mgt->Add(graph_data_altTag);
        mgt->Draw("ap");
        mgt->GetXaxis()->SetTitle("relIso");
        mgt->GetYaxis()->SetTitle("Efficiency "+namesPt+", "+gribu);
        legendt->AddEntry(graph_data_altTag,"data");
        legendt->AddEntry(graph_MC_altTag,"MC");
        legendt->Draw();
        oit->SaveAs("plots/"+folderName+"/UL_eta_"+gribu+"_"+namesPt+"_Efficiency_"+gribuTag+"_altTag.png");
      
      }
    }
  }
}

