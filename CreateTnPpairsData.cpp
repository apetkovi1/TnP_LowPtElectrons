#include "TTree.h"
#include "TChain.h"
#include "TFile.h"

void CreateTnPpairsData()
{
  float TagPt,ProbePt,Diele_mass,TagEta,ProbeEta;
  bool ProbePass;
  std::vector<float> *ele_pt=0;
  std::vector<float> *ele_eta=0;
  std::vector<bool> *mvaEleID_Fall17_noIso_V2_wp90=0, *mvaEleID_Fall17_noIso_V2_wp80=0;
  TFile* file = TFile::Open("Bpark_DATA_2018.root");
  TTree* originalTree = (TTree*)file->Get("electrons/Events");
  originalTree->SetBranchStatus("*",0);
  originalTree->SetBranchStatus("Diele_mass",1);
  originalTree->SetBranchStatus("ele_pt",1);
  originalTree->SetBranchStatus("mvaEleID_Fall17_noIso_V2_wp90",1);
  originalTree->SetBranchStatus("mvaEleID_Fall17_noIso_V2_wp80",1);
  originalTree->SetBranchStatus("ele_eta",1);
  TFile* output = TFile::Open("TnPpairs_DATA.root","RECREATE");
  TTree* selectedTree = originalTree->CopyTree("mvaEleID_Fall17_noIso_V2_wp90[0]==1 || mvaEleID_Fall17_noIso_V2_wp90[1]==1");
  selectedTree->SetBranchAddress("Diele_mass",&Diele_mass);
  selectedTree->SetBranchAddress("ele_pt",&ele_pt);
  selectedTree->SetBranchAddress("mvaEleID_Fall17_noIso_V2_wp90",&mvaEleID_Fall17_noIso_V2_wp90);
  selectedTree->SetBranchAddress("mvaEleID_Fall17_noIso_V2_wp80",&mvaEleID_Fall17_noIso_V2_wp80);
  selectedTree->SetBranchAddress("ele_eta",&ele_eta);
  TBranch *BranchTagPt=selectedTree->Branch("TagPt",&TagPt);
  TBranch *BranchProbePt=selectedTree->Branch("ProbePt",&ProbePt);
  TBranch *BranchTagEta=selectedTree->Branch("TagEta",&TagEta);
  TBranch *BranchProbeEta=selectedTree->Branch("ProbeEta",&ProbeEta);
  TBranch *BranchProbePass=selectedTree->Branch("ProbePass",&ProbePass);
  Long64_t nentries = selectedTree->GetEntries();
  for(Long64_t i=0;i<nentries;i++)
  {
    selectedTree->GetEntry(i);
    if(mvaEleID_Fall17_noIso_V2_wp90->at(0)==1 && mvaEleID_Fall17_noIso_V2_wp80->at(1)==1)
    {
      TagPt=ele_pt->at(0);
      ProbePt=ele_pt->at(1);
      TagEta=ele_eta->at(0);
      ProbeEta=ele_eta->at(1);
      ProbePass=1;
    }
    if(mvaEleID_Fall17_noIso_V2_wp90->at(0)==1 && mvaEleID_Fall17_noIso_V2_wp80->at(1)==0)
    {
      TagPt=ele_pt->at(0);
      ProbePt=ele_pt->at(1);
      TagEta=ele_eta->at(0);
      ProbeEta=ele_eta->at(1);
      ProbePass=0;
    }
    if(mvaEleID_Fall17_noIso_V2_wp80->at(0)==0 && mvaEleID_Fall17_noIso_V2_wp90->at(1)==1)
    {
      TagPt=ele_pt->at(1);
      ProbePt=ele_pt->at(0);
      TagEta=ele_eta->at(1);
      ProbeEta=ele_eta->at(0);
      ProbePass=0;
    }

    BranchTagPt->Fill();
    BranchTagEta->Fill();
    BranchProbeEta->Fill();
    BranchProbePt->Fill();
    BranchProbePass->Fill();
  }
  selectedTree->Write();

}
