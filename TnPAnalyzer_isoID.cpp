#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "RooCrystalBall.cxx"
#include "RooCrystalBall.h"
#include "MakeFitData.cpp"
#include "MakeFitDataAltBgk.cpp"
#include "MakeFitDataAltSig.cpp"
#include "MakeFitDataAltSigBkg.cpp"
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

void TnPAnalyzer_isoID()
{
  // Some definitions
  int i;
  float ptBins[]={5,7,10,15};  //this part is configurable, select bins for fitting
  string tag_cut_status="tagNoIso"; //decoration for png names
  string folderName="iso_HZZ"; //ID Choice;
  bool save_status=1;
  bool doEfficiency=1;
  int variations = 4;
  bool FitInnerBarrel=1, FitOuterBarrel=0, FitEndcap=0;
  TString EtaCut="";
  TString isoCut="";
  TString ptCut="";
  int etaBin=0;
  int isoBin=0;
  int pTbin=0;
  TString etaBinName="";

  string outFileName="plots/"+folderName+"/UL_altRez_eff.txt"; // Efficiency output file name
  // creating output file to store the results
  std::ofstream outfile;
  outfile.open(outFileName, std::ios_base::app); // append instead of overwrite
  
  outfile<<"eta\tpt\trelIso\td_eff_n\td_eff_n_err\tm_eff\tm_eff_err\td_eff_b\td_eff_b_err\td_eff_s\td_eff_s_err\td_eff_sb\td_eff_sb_err\td_eff_mean\td_eff_mean_err\tcomment"<<std::endl;  // Column names in output file
  //outfile<<"eta\tpt\trelIso\td_eff_n\td_eff_n_err\tm_eff\tm_eff_err\td_eff_b\td_eff_b_err\td_eff_s\td_eff_s_err\td_eff_t\td_eff_t_err\tm_eff_t\tm_eff_t_err\tcomment"<<std::endl; // altTag case commneted out
  float MassConstrainPass, MassConstrainFail;

  //float mvaCutValue=0.97; //altTag   // altTag case commneted out
  //TString altCutData=TString::Format(" TagMVA>%f ",mvaCutValue); // altTag case commneted out
  //TString altCutMC=TString::Format(" TagMVA>%f && TagPt>7 ",mvaCutValue); // altTag case commneted out
  
  
  
  
  // importing DATA file
  TFile *file_data = TFile::Open("/eos/user/a/angaile/bParkingRootFiles/wpHZZSum18/TnPpairs_UL_Data_wpHZZSum.root"); // Input TnP_ntuple(Output from "CreateTnPpairsData.cpp") for DATA and define which branches to use
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
  TFile *file_MC = TFile::Open("/eos/user/a/angaile/bParkingRootFiles/wpHZZSum18/TnPpairs_UL_MC_wpHZZSum.root"); // Input TnP_ntuple(Output from "CreateTnPpairsMC.cpp") for MC and define which branches to use
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
  
  for(etaBin=1;etaBin<=3;etaBin++) // iterate through eta bins
	{
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
    
    std::cout <<EtaCut << std::endl; 
    int NumberOfPtBins=sizeof(ptBins)/sizeof(ptBins[0])-1;
    
    for(pTbin=0;pTbin<NumberOfPtBins;pTbin++)
    {  
      ptCut=TString::Format(" ProbePt>%f && ProbePt<%f ",ptBins[pTbin],ptBins[pTbin+1]);
      // for each eta region
      // Define (eff, eff_error) pairs for all fit models. (MC, DATA -> Nominal, altSig, altBkg, altSigBkg)
      std::vector<std::pair<float,float>> YieldErrorDataPass, YieldErrorDataFail, YieldErrorMCPass, YieldErrorMCFail;
      std::vector<std::pair<float,float>> YieldErrorDataPassBkg, YieldErrorDataFailBkg, YieldErrorDataPassSig, YieldErrorDataFailSig, YieldErrorDataPassSigBkg, YieldErrorDataFailSigBkg;
      //std::vector<std::pair<float,float>> YieldErrorDataPassTag, YieldErrorDataFailTag, YieldErrorMCPassTag, YieldErrorMCFailTag; // altTag case commneted out

      // Defining electron's relative isolation bins in which the T&P efficiency measurements will be done.
      //Bin boundaries are defined by the analyzer himself/herself, based on the rel_Iso distribution study. (Define bins in a way that there are similar amount of events in each [pT/eta] bin). For example, here the bins are defined in a way that there are similar amount of failling probes in each bin, for the iso_HZZ ID.
      float isoBins[5]={0};
      if(ptBins[pTbin+1]== 7 )
      {
        if(etaBin==1) {
          float pisoBins[]={0,0.001,0.155,0.330,0.5};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
        if(etaBin==2) {
          float pisoBins[]={0,0.001,0.155,0.324,0.5};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
        if(etaBin==3) {
          float pisoBins[]={0,0.001,0.204,0.5,2};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
      } else if(ptBins[pTbin+1]== 10 )
      {
        if(etaBin==1) {
          float pisoBins[]={0,0.081,0.277,0.452,0.6};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
        if(etaBin==2) {
          float pisoBins[]={0,0.095,0.286,0.455,0.6};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
        if(etaBin==3) {
          float pisoBins[]={0,0.065,0.293,0.6,2};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
      } else if(ptBins[pTbin+1]== 15 )
      {
        if(etaBin==1) {
          float pisoBins[]={0,0.040,0.141,0.226,0.3};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
        if(etaBin==2) {
          float pisoBins[]={0,0.037,0.139,0.224,0.3};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
        if(etaBin==3) {
          float pisoBins[]={0,0.040,0.164,0.3,2};
          memcpy(isoBins, pisoBins, sizeof(isoBins));
        }
      }

      int NumberOfIsoBins=sizeof(isoBins)/sizeof(isoBins[0])-1;
      int NumberOfBins=NumberOfIsoBins; //NumberOfPtBins*NumberOfIsoBins;

      // Passing and failling probes (data and MC) // Define trees in which binned data will be held for fits
      TTree* PassingProbesData[NumberOfBins];
      TTree* FaillingProbesData[NumberOfBins];
      TTree* PassingProbesMC[NumberOfBins];
      TTree* FaillingProbesMC[NumberOfBins];
      
      // Passing and failling probes for alt Tag (data and MC)
      //TTree* tPassingProbesData[NumberOfBins]; // altTag case commneted out
      //TTree* tFaillingProbesData[NumberOfBins]; // altTag case commneted out
      //TTree* tPassingProbesMC[NumberOfBins]; // altTag case commneted out
      //TTree* tFaillingProbesMC[NumberOfBins]; // altTag case commneted out
      
      string namesBins[NumberOfBins];
      string namesPt="pt_"+std::to_string((int)ptBins[pTbin])+"_"+std::to_string((int)ptBins[pTbin+1]);
      i=0;
      for(isoBin=0;isoBin<NumberOfIsoBins;isoBin++) // Iterating through rel_iso bins // // Fill the trees with binned data
      {
        isoCut=TString::Format(" ProbeRelISO>=%f && ProbeRelISO<%f ",isoBins[isoBin],isoBins[isoBin+1]);
      
        namesBins[i] = namesPt+"_relIso_"+std::to_string((int)(isoBins[isoBin]*100))+"_"+std::to_string((int)(isoBins[isoBin+1]*100));
      
        PassingProbesData[i] = tree_data->CopyTree(ptCut+"&&"+isoCut + " && ProbePass==1 && " + EtaCut);
        FaillingProbesData[i] = tree_data->CopyTree(ptCut+"&&"+isoCut + " && ProbePass==0 && " + EtaCut);
        
        PassingProbesMC[i] = tree_MC->CopyTree(ptCut+"&&"+isoCut + " && ProbePass==1 && " + EtaCut + " && TagPt>7");
        FaillingProbesMC[i] = tree_MC->CopyTree(ptCut+"&&"+isoCut + " && ProbePass==0 && " + EtaCut + " && TagPt>7");
         
        //tPassingProbesData[i] = tree_data->CopyTree(ptCut+"&&"+isoCut + " && ProbePass==1 && " + EtaCut+" && "+ altCutData);  // altTag case commneted out
        //tFaillingProbesData[i] = tree_data->CopyTree(ptCut+"&&"+isoCut + " && ProbePass==0 && " + EtaCut+" && "+ altCutData); // altTag case commneted out
        //tPassingProbesMC[i] = tree_MC->CopyTree(ptCut+"&&"+isoCut + " && ProbePass==1 && " + EtaCut+" && "+ altCutMC); // altTag case commneted out
        //tFaillingProbesMC[i] = tree_MC->CopyTree(ptCut+"&&"+isoCut + " && ProbePass==0 && "+ EtaCut+" && "+ altCutMC); // altTag case commneted out
        i++;
      }

      //Fit for passing and failling probes in each bin
      for(i=0;i<NumberOfBins;i++) // 
      {
        std::pair<float,float> InfoMCPass, InfoMCFail; // , InfoMCPassT, InfoMCFailT; // altTag case commneted out
        std::tuple<float,float,float,float,float,float,float,float,float> ParDataPass, ParDataFail; // ,  ParDataPassT, ParDataFailT; // altTag case commneted out

        // Get info about DSCB parameters from MC fit -> Will be used for DATA fit as initial values and boundaries
        InfoMCPass=GetFitMC(PassingProbesMC[i]);
        ParDataPass=GetFitData(PassingProbesData[i],InfoMCPass);
        InfoMCFail=GetFitMC(FaillingProbesMC[i]);
        ParDataFail=GetFitData(FaillingProbesData[i],InfoMCFail);

        //InfoMCPassT=GetFitMC(tPassingProbesMC[i]); // altTag case commneted out
        //ParDataPassT=GetFitData(tPassingProbesData[i],InfoMCPassT); // altTag case commneted out
        //InfoMCFailT=GetFitMC(tFaillingProbesMC[i]); // altTag case commneted out
        //ParDataFailT=GetFitData(tFaillingProbesData[i],InfoMCFailT); // altTag case commneted out

      
        string plotNamePart="plots/"+folderName+"/UL_eta_"+std::to_string(etaBin)+"_"+namesBins[i]+"_"+tag_cut_status;
        if(isoBins[i+1]==2) // Do not save 4th bin in endcap, as there are only 3 valid bins, due to limited singla statistics. 4th bin boundary defined as value 2 to identify it here.
          save_status=0;
        else
          save_status=1;
        // regular fits
        YieldErrorDataPass.push_back(MakeFitDataNom(PassingProbesData[i],ParDataPass, save_status, plotNamePart+"_nom_PassingProbesData.png",false));
        //std::cout << "it was passing data " << namesBins[i] << "\n\n\n\n\n";
        YieldErrorDataFail.push_back(MakeFitDataNom(FaillingProbesData[i],ParDataFail, save_status, plotNamePart+"_nom_FaillingProbesData.png",false));
        //std::cout << "it was failing data " << namesBins[i] << "\n\n\n\n\n";    
        YieldErrorMCPass.push_back(MakeFitMC(PassingProbesMC[i], save_status, plotNamePart+"_nom_PassingProbesMC.png"));
        YieldErrorMCFail.push_back(MakeFitMC(FaillingProbesMC[i], save_status, plotNamePart+"_nom_FaillingProbesMC.png"));
      
        // alt tag fits // altTag case commneted out
        //YieldErrorDataPassTag.push_back(MakeFitDataNom(tPassingProbesData[i],ParDataPassT, save_status,  plotNamePart+"_altTag_PassingProbesData.png",false)); // altTag case commneted out
        //YieldErrorDataFailTag.push_back(MakeFitDataNom(tFaillingProbesData[i],ParDataFailT, save_status, plotNamePart+"_altTag_FaillingProbesData.png",false)); // altTag case commneted out
        //YieldErrorMCPassTag.push_back(MakeFitMC(tPassingProbesMC[i], save_status, plotNamePart+"_altTag_PassingProbesMC.png")); // altTag case commneted out
        //YieldErrorMCFailTag.push_back(MakeFitMC(tFaillingProbesMC[i], save_status, plotNamePart+"_altTag_FaillingProbesMC.png")); // altTag case commneted out
      
        // alt bgk fits
        YieldErrorDataPassBkg.push_back(MakeFitDataAltBgk(PassingProbesData[i],ParDataPass, save_status,  plotNamePart+"_altBgk_PassingProbesData.png",false));
        YieldErrorDataFailBkg.push_back(MakeFitDataAltBgk(FaillingProbesData[i],ParDataFail, save_status, plotNamePart+"_altBgk_FaillingProbesData.png",false));
        
        // alt sign fits
        //std::cout << "\n\n\nStarting alt\n\n\n "  << std::endl;
        YieldErrorDataPassSig.push_back(MakeFitDataAltSig(PassingProbesData[i],ParDataPass, save_status,  plotNamePart+"_altSig_PassingProbesData.png",false));
        YieldErrorDataFailSig.push_back(MakeFitDataAltSig(FaillingProbesData[i],ParDataFail, save_status, plotNamePart+"_altSig_FaillingProbesData.png",false));

        // alt sig alt bkg fits
        //std::cout << "\n\n\nStarting alt\n\n\n "  << std::endl;
        YieldErrorDataPassSigBkg.push_back(MakeFitDataAltSigBkg(PassingProbesData[i],ParDataPass, save_status,  plotNamePart+"_altSigBkg_PassingProbesData.png", false));
        YieldErrorDataFailSigBkg.push_back(MakeFitDataAltSigBkg(FaillingProbesData[i],ParDataFail, save_status, plotNamePart+"_altSigBkg_FaillingProbesData.png", false));
      

      }

      // here it should just be relIso
      TEfficiency* Eff_MC = Efficiency(YieldErrorMCPass, YieldErrorMCFail, isoBins, NumberOfBins);
      //TEfficiency* Eff_MC_altTag = Efficiency(YieldErrorMCPassTag, YieldErrorMCFailTag, isoBins, NumberOfBins); // altTag case commneted out

      TEfficiency* Eff_DATA = Efficiency(YieldErrorDataPass, YieldErrorDataFail, isoBins, NumberOfBins);
      // std::cout << "it was Eff calc " << "\n";  
      TEfficiency* Eff_DATA_altBkg = Efficiency(YieldErrorDataPassBkg, YieldErrorDataFailBkg, isoBins, NumberOfBins);
      TEfficiency* Eff_DATA_altSig = Efficiency(YieldErrorDataPassSig, YieldErrorDataFailSig, isoBins, NumberOfBins);
      //TEfficiency* Eff_DATA_altTag = Efficiency(YieldErrorDataPassTag, YieldErrorDataFailTag, isoBins, NumberOfBins); // altTag case commneted out
      TEfficiency* Eff_DATA_altSigBkg = Efficiency(YieldErrorDataPassSigBkg, YieldErrorDataFailSigBkg, isoBins, NumberOfBins);

      


      
      for(i=1;i<=NumberOfIsoBins;i++) // Write efficiency values into output file.
      {
        if(isoBins[i+1]==2)
          continue;

        double Eff_DATA_mean = (Eff_DATA->GetEfficiency(i)+Eff_DATA_altBkg->GetEfficiency(i)+Eff_DATA_altSig->GetEfficiency(i)+Eff_DATA_altSigBkg->GetEfficiency(i))/4;
        double Eff_DATA_RMS = sqrt(pow(Eff_DATA->GetEfficiency(i)-Eff_DATA_mean,2)+pow(Eff_DATA_altBkg->GetEfficiency(i)-Eff_DATA_mean,2)+pow(Eff_DATA_altSig->GetEfficiency(i)-Eff_DATA_mean,2)+pow(Eff_DATA_altSigBkg->GetEfficiency(i)-Eff_DATA_mean,2)/(variations-1));
        double Eff_DATA_mean_err = sqrt(pow(Eff_DATA_RMS/sqrt(variations),2)+pow(Eff_DATA->GetEfficiencyErrorUp(i),2));

        outfile<<std::to_string(etaBin)<<"\t"<<std::to_string((int)ptBins[pTbin+1])<<"\t"<<std::to_string((float)isoBins[i-1])<<"\t"<<Eff_DATA->GetEfficiency(i)<<"\t"<<Eff_DATA->GetEfficiencyErrorUp(i)<<"\t"<<Eff_MC->GetEfficiency(i)<<"\t"<<Eff_MC->GetEfficiencyErrorUp(i)<<"\t"<<Eff_DATA_altBkg->GetEfficiency(i)<<"\t"<<Eff_DATA_altBkg->GetEfficiencyErrorUp(i)<<"\t"<<Eff_DATA_altSig->GetEfficiency(i)<<"\t"<<Eff_DATA_altSig->GetEfficiencyErrorUp(i)<<"\t"<<Eff_DATA_altSigBkg->GetEfficiency(i)<<"\t"<<Eff_DATA_altSigBkg->GetEfficiencyErrorUp(i)<<"\t"<<Eff_DATA_mean<<"\t"<<Eff_DATA_mean_err<<"\t"<<namesBins[i-1]<<std::endl;
      }
      
      //Save eff object to root file
      //  TString effName="plots/"+idChoice+"/UL_eta_"+etaBinName+"_"+namesPt+"_efficiency_"+tag_cut_status+".root";
      TFile* eff = TFile::Open("plots/"+folderName+"/UL_eta_"+etaBinName+"_"+namesPt+"_efficiency_"+tag_cut_status+".root","RECREATE");
      eff->WriteObject(Eff_MC, "Eff_MC");
      //eff->WriteObject(Eff_MC_altTag, "Eff_MC_altTag"); // altTag case commneted out
      eff->WriteObject(Eff_DATA, "Eff_DATA");
      eff->WriteObject(Eff_DATA_altSig, "Eff_DATA_altSig");
      eff->WriteObject(Eff_DATA_altBkg, "Eff_DATA_altBkg");
      //eff->WriteObject(Eff_DATA_altTag, "Eff_DATA_altTag"); // altTag case commneted out
      eff->WriteObject(Eff_DATA_altSigBkg, "Eff_DATA_altSigBkg");
      

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
        mg->GetYaxis()->SetTitle("Efficiency "+namesPt+", "+etaBinName);
        legend->AddEntry(graph_data,"data");
        legend->AddEntry(graph_MC,"MC");
        legend->Draw();
        oi->SaveAs("plots/"+folderName+"/UL_eta_"+etaBinName+"_"+namesPt+"_Efficiency_"+tag_cut_status+"_nom.png");

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
        mgb->GetYaxis()->SetTitle("Efficiency_altBkg "+namesPt+", "+etaBinName);
        legendb->AddEntry(graph_datab,"data");
        legendb->AddEntry(graph_MCb,"MC");
        legendb->Draw();
        oib->SaveAs("plots/"+folderName+"/UL_eta_"+etaBinName+"_"+namesPt+"_Efficiency_"+tag_cut_status+"_altBkg.png");
        
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
        mgs->GetYaxis()->SetTitle("Efficiency_altSig "+namesPt+", "+etaBinName);
        legends->AddEntry(graph_datas,"data");
        legends->AddEntry(graph_MCs,"MC");
        legends->Draw();
        ois->SaveAs("plots/"+folderName+"/UL_eta_"+etaBinName+"_"+namesPt+"_Efficiency_"+tag_cut_status+"_altSig.png");
        
        //efficiency plot alternative SigBkg
        auto legendt = new TLegend(0.9,0.9,1.0,1.0);
        auto mgt  = new TMultiGraph();
        Eff_MC->SetLineColor(kRed);
        Eff_DATA_altSigBkg->SetLineColor(kBlue);
        TCanvas* oit = new TCanvas();
        oit->cd();
        Eff_MC->Draw();
        gPad->Update();
        auto graph_MCsb = Eff_MC->GetPaintedGraph();
        mgt->Add(graph_MCsb);
        Eff_DATA_altSigBkg->Draw();
        gPad->Update();
        auto graph_data_altSigBkg = Eff_DATA_altSigBkg->GetPaintedGraph();
        mgt->Add(graph_data_altSigBkg);
        mgt->Draw("ap");
        mgt->GetXaxis()->SetTitle("relIso");
        mgt->GetYaxis()->SetTitle("Efficiency "+namesPt+", "+etaBinName);
        legendt->AddEntry(graph_data_altSigBkg,"data");
        legendt->AddEntry(graph_MCsb,"MC");
        legendt->Draw();
        oit->SaveAs("plots/"+folderName+"/UL_eta_"+etaBinName+"_"+namesPt+"_Efficiency_"+tag_cut_status+"_altSigBkg.png");
      
      }
    }
  }
}

