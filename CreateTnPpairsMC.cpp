#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include <stdlib.h>

void CreateTnPpairsMC()
{
  int j,TagIndex,ProbeIndex;
  float TagPt,ProbePt,Diele_mass,TagEta,ProbeEta, dR;
  bool ProbePass;
  std::vector<float> *ele_pt=0, *ele_eta=0,*ele_phi=0, *ElectronMVAEstimatorRun2Fall17NoIsoV2Values=0, *ElectronMVAEstimatorRun2Fall17IsoV2Values=0;
  std::vector<bool> *mvaEleID_Fall17_noIso_V2_wp90=0, *mvaEleID_Fall17_noIso_V2_wp80=0, *mvaEleID_Fall17_noIso_V2_wpLoose_unsopported=0,
  *mvaEleID_Fall17_iso_V2_wpHZZ_unsopported=0;
  std::vector<int> *ele_mother=0;

  srand (time(NULL));

  TFile* file = TFile::Open("Bpark_MC_v2.root");
  TTree* originalTree = (TTree*)file->Get("electrons/Events");
  originalTree->SetBranchStatus("*",0);
  originalTree->SetBranchStatus("Diele_mass",1);
  originalTree->SetBranchStatus("ele_pt",1);
  originalTree->SetBranchStatus("mvaEleID_Fall17_noIso_V2_wp90",1);
  originalTree->SetBranchStatus("mvaEleID_Fall17_noIso_V2_wp80",1);
  originalTree->SetBranchStatus("mvaEleID_Fall17_noIso_V2_wpLoose_unsopported",1);
  originalTree->SetBranchStatus("mvaEleID_Fall17_iso_V2_wpHZZ_unsopported",1);
  originalTree->SetBranchStatus("ele_eta",1);
  originalTree->SetBranchStatus("ElectronMVAEstimatorRun2Fall17NoIsoV2Values",1);
  originalTree->SetBranchStatus("ElectronMVAEstimatorRun2Fall17IsoV2Values",1);
  originalTree->SetBranchStatus("ele_mother",1);
  originalTree->SetBranchStatus("ele_phi",1);
  originalTree->SetBranchStatus("dR",1);
  TFile* output = TFile::Open("TnPpairs_MC.root","RECREATE");

  TTree* selectedTree = originalTree->CopyTree("ele_mother[0]==443 && ele_mother[1]==443 && (ElectronMVAEstimatorRun2Fall17NoIsoV2Values[0]>0.95 && ele_pt[0]>7 && (fabs(ele_eta[0])<1.44 || fabs(ele_eta[0])>1.57)) || (ElectronMVAEstimatorRun2Fall17NoIsoV2Values[1]>0.95 && ele_pt[1]>7 && (fabs(ele_eta[1])<1.44 || fabs(ele_eta[1])>1.57))");
  selectedTree->SetBranchAddress("Diele_mass",&Diele_mass);
  selectedTree->SetBranchAddress("ele_pt",&ele_pt);
  selectedTree->SetBranchAddress("mvaEleID_Fall17_noIso_V2_wp90",&mvaEleID_Fall17_noIso_V2_wp90);
  selectedTree->SetBranchAddress("mvaEleID_Fall17_noIso_V2_wp80",&mvaEleID_Fall17_noIso_V2_wp80);
  selectedTree->SetBranchAddress("mvaEleID_Fall17_noIso_V2_wpLoose_unsopported",&mvaEleID_Fall17_noIso_V2_wpLoose_unsopported);
  selectedTree->SetBranchAddress("mvaEleID_Fall17_iso_V2_wpHZZ_unsopported",&mvaEleID_Fall17_iso_V2_wpHZZ_unsopported);
  selectedTree->SetBranchAddress("ele_eta",&ele_eta);
  selectedTree->SetBranchAddress("ElectronMVAEstimatorRun2Fall17NoIsoV2Values",&ElectronMVAEstimatorRun2Fall17NoIsoV2Values);
  selectedTree->SetBranchAddress("ElectronMVAEstimatorRun2Fall17IsoV2Values",&ElectronMVAEstimatorRun2Fall17IsoV2Values);
  selectedTree->SetBranchAddress("ele_mother",&ele_mother);
  selectedTree->SetBranchAddress("ele_phi",&ele_phi);
  selectedTree->SetBranchAddress("dR",&dR);
  TBranch *BranchTagPt=selectedTree->Branch("TagPt",&TagPt);
  TBranch *BranchProbePt=selectedTree->Branch("ProbePt",&ProbePt);
  TBranch *BranchTagEta=selectedTree->Branch("TagEta",&TagEta);
  TBranch *BranchProbeEta=selectedTree->Branch("ProbeEta",&ProbeEta);
  TBranch *BranchProbePass=selectedTree->Branch("ProbePass",&ProbePass);
  TBranch *BranchProbeIndex=selectedTree->Branch("ProbeIndex",&ProbeIndex);
  TBranch *BranchTagIndex=selectedTree->Branch("TagIndex",&TagIndex);
  Long64_t nentries = selectedTree->GetEntries();
  for(Long64_t i=0;i<nentries;i++) //Tag-Probe pairs, Probe-Tag pairs and Tag-Tag pairs in which first is selected as tag
  {
    selectedTree->GetEntry(i);
    if((ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(0)>0.95 && ele_pt->at(0)>7 && (fabs(ele_eta->at(0))<1.44 || fabs(ele_eta->at(0))>1.57)) && (ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(1)>0.95 && ele_pt->at(1)>7 &&(fabs(ele_eta->at(1))<1.44 || fabs(ele_eta->at(1))>1.57)))
    TagIndex=0;
    else
    for(j=0;j<2;j++)
    {
      if(ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(j)>0.95 && ele_pt->at(j)>7 && (fabs(ele_eta->at(j))<1.44 || fabs(ele_eta->at(j))>1.57)) //Tag selection
      TagIndex=j;
    }
    ProbeIndex=abs(TagIndex-1);

    TagPt=ele_pt->at(TagIndex);
    ProbePt=ele_pt->at(ProbeIndex);
    TagEta=ele_eta->at(TagIndex);
    ProbeEta=ele_eta->at(ProbeIndex);
    if(mvaEleID_Fall17_noIso_V2_wpLoose_unsopported->at(ProbeIndex)==1) //ID to measure efficiency
    ProbePass=1;
    else
    ProbePass=0;
    BranchTagPt->Fill();
    BranchTagEta->Fill();
    BranchProbeEta->Fill();
    BranchProbePt->Fill();
    BranchProbePass->Fill();
    BranchTagIndex->Fill();
    BranchProbeIndex->Fill();
  }

  TTree* TempTree = originalTree->CopyTree("ele_mother[0]==443 && ele_mother[1]==443 && (ElectronMVAEstimatorRun2Fall17NoIsoV2Values[0]>0.95 && ele_pt[0]>7 && (fabs(ele_eta[0])<1.44 || fabs(ele_eta[0])>1.57)) && (ElectronMVAEstimatorRun2Fall17NoIsoV2Values[1]>0.95 && ele_pt[1]>7 && (fabs(ele_eta[1])<1.44 || fabs(ele_eta[1])>1.57))");
  TempTree->SetBranchAddress("Diele_mass",&Diele_mass);
  TempTree->SetBranchAddress("ele_pt",&ele_pt);
  TempTree->SetBranchAddress("mvaEleID_Fall17_noIso_V2_wp90",&mvaEleID_Fall17_noIso_V2_wp90);
  TempTree->SetBranchAddress("mvaEleID_Fall17_noIso_V2_wp80",&mvaEleID_Fall17_noIso_V2_wp80);
  TempTree->SetBranchAddress("mvaEleID_Fall17_noIso_V2_wpLoose_unsopported",&mvaEleID_Fall17_noIso_V2_wpLoose_unsopported);
  TempTree->SetBranchAddress("mvaEleID_Fall17_iso_V2_wpHZZ_unsopported",&mvaEleID_Fall17_iso_V2_wpHZZ_unsopported);
  TempTree->SetBranchAddress("ele_eta",&ele_eta);
  TempTree->SetBranchAddress("ElectronMVAEstimatorRun2Fall17NoIsoV2Values",&ElectronMVAEstimatorRun2Fall17NoIsoV2Values);
  TempTree->SetBranchAddress("ElectronMVAEstimatorRun2Fall17IsoV2Values",&ElectronMVAEstimatorRun2Fall17IsoV2Values);
  TempTree->SetBranchAddress("ele_mother",&ele_mother);
  TBranch *BranchTagPt_=TempTree->Branch("TagPt",&TagPt);
  TBranch *BranchProbePt_=TempTree->Branch("ProbePt",&ProbePt);
  TBranch *BranchTagEta_=TempTree->Branch("TagEta",&TagEta);
  TBranch *BranchProbeEta_=TempTree->Branch("ProbeEta",&ProbeEta);
  TBranch *BranchProbePass_=TempTree->Branch("ProbePass",&ProbePass);
  TBranch *BranchProbeIndex_=TempTree->Branch("ProbeIndex",&ProbeIndex);
  TBranch *BranchTagIndex_=TempTree->Branch("TagIndex",&TagIndex);
  nentries = TempTree->GetEntries();
  for(Long64_t i=0;i<nentries;i++) //Tag-Tag pairs in which second is selected as tag
  {
    TempTree->GetEntry(i);
    TagIndex=1;
    ProbeIndex=0;
    TagPt=ele_pt->at(TagIndex);
    ProbePt=ele_pt->at(ProbeIndex);
    TagEta=ele_eta->at(TagIndex);
    ProbeEta=ele_eta->at(ProbeIndex);
    if(mvaEleID_Fall17_noIso_V2_wpLoose_unsopported->at(ProbeIndex)==1) //ID to measure efficiency
    ProbePass=1;
    else
    ProbePass=0;
    BranchTagPt_->Fill();
    BranchTagEta_->Fill();
    BranchProbeEta_->Fill();
    BranchProbePt_->Fill();
    BranchProbePass_->Fill();
    BranchTagIndex_->Fill();
    BranchProbeIndex_->Fill();
  }

  TList *list = new TList;
  list->Add(selectedTree);
  list->Add(TempTree);
  TTree *FinalTree = TTree::MergeTrees(list);
  FinalTree->Write();
}
