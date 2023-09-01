#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include <stdlib.h>

void CreateTnPpairsMC()
{
  // Define parameters and branches used (Multiple lines needs to added, check the already defined branches and define new in the same manner) (Example "ele_pt" and corresponding TagPt/ProbePt branches/parameters)
  // If some of the defined branches here are not included in ntuples, comment them out of this code
  int j,TagIndex,ProbeIndex,ele_sameVertex,ele_ip3D_match;
  float TagPt,ProbePt,Diele_mass,Diele_pt,TagEta,ProbeEta,TagPhi,ProbePhi,dR,TagRelISO,ProbeRelISO,ProbeChIso,ProbeNeuIso,ProbePhoIso,TagDxy,ProbeDxy,TagDxyError,ProbeDxyError,TagDz,ProbeDz,TagDzError,ProbeDzError,TagMVA;
  bool ProbePass;
  std::vector<float> *ele_pt=0, *ele_eta=0, *ele_dxy=0, *ele_dz=0, *ele_dxyError=0, *ele_dzError=0, *ele_phi=0, *ele_pfPhotonIso=0 , *ele_pfChargedHadIso=0 , *ele_pfNeutralHadIso=0 , *ElectronMVAEstimatorRun2Fall17NoIsoV2Values=0, *ElectronMVAEstimatorRun2Fall17IsoV2Values=0,*relISO_a_corr=0;
  std::vector<bool> *mvaEleID_Fall17_noIso_V2_wp90=0, *mvaEleID_Fall17_noIso_V2_wp80=0, *mvaEleID_Fall17_noIso_V2_wpLoose_unsopported=0,
  *mvaEleID_Fall17_iso_V2_wpHZZ_unsopported=0,*mvaEleID_Fall17_iso_V2_wp80=0,*mvaEleID_Fall17_iso_V2_wp90=0;
  std::vector<int> *ele_mother=0;

  srand (time(NULL));


  TFile* file = TFile::Open("/eos/user/n/nstrautn/Bparking_DATA/Bpark_22_SM_test_MC.root"); // Change input MC ntuple file here (Ntuple should come from the related miniAOD Ntuplizer)
  TTree* originalTree = (TTree*)file->Get("electrons/Events");
  originalTree->SetBranchStatus("*",0);
  originalTree->SetBranchStatus("Diele_mass",1);
  originalTree->SetBranchStatus("Diele_pt",1);
  originalTree->SetBranchStatus("dR",1);
  originalTree->SetBranchStatus("relISO_a_corr",1);
  originalTree->SetBranchStatus("ele_pt",1);
  originalTree->SetBranchStatus("ele_sameVertex",1);
  originalTree->SetBranchStatus("ele_ip3D_match",1);
  originalTree->SetBranchStatus("ele_dxy",1);
  originalTree->SetBranchStatus("ele_dxyError",1);
  originalTree->SetBranchStatus("ele_dz",1);
  originalTree->SetBranchStatus("ele_dzError",1);
  originalTree->SetBranchStatus("ele_phi",1);
  originalTree->SetBranchStatus("ele_pfPhotonIso",1);
  originalTree->SetBranchStatus("ele_pfChargedHadIso",1);
  originalTree->SetBranchStatus("ele_pfNeutralHadIso",1);
  originalTree->SetBranchStatus("mvaEleID_Fall17_noIso_V2_wp90",1);
  originalTree->SetBranchStatus("mvaEleID_Fall17_noIso_V2_wp80",1);
  originalTree->SetBranchStatus("mvaEleID_Fall17_iso_V2_wp90",1);
  originalTree->SetBranchStatus("mvaEleID_Fall17_iso_V2_wp80",1);
  originalTree->SetBranchStatus("mvaEleID_Fall17_noIso_V2_wpLoose_unsopported",1);
  originalTree->SetBranchStatus("mvaEleID_Fall17_iso_V2_wpHZZ_unsopported",1);
  originalTree->SetBranchStatus("ele_eta",1);
  originalTree->SetBranchStatus("ElectronMVAEstimatorRun2Fall17NoIsoV2Values",1);
  originalTree->SetBranchStatus("ElectronMVAEstimatorRun2Fall17IsoV2Values",1);
  originalTree->SetBranchStatus("ele_mother",1);
  TFile* output = TFile::Open("TnPpairs_MC.root","RECREATE");

  // Define tag and probe tag cuts and pair selection (Check all places where this cut is written. Change needed in 4 lines) 1/4
  TTree* selectedTree = originalTree->CopyTree("ele_mother[0]==443 && ele_mother[1]==443 && (ElectronMVAEstimatorRun2Fall17NoIsoV2Values[0]>0.95 && ele_pt[0]>7 && (fabs(ele_eta[0])<1.44 || fabs(ele_eta[0])>1.57)) || (ElectronMVAEstimatorRun2Fall17NoIsoV2Values[1]>0.95 && ele_pt[1]>7 && (fabs(ele_eta[1])<1.44 || fabs(ele_eta[1])>1.57))"); 
  selectedTree->SetBranchAddress("Diele_mass",&Diele_mass);
  selectedTree->SetBranchAddress("Diele_pt",&Diele_pt);
  selectedTree->SetBranchAddress("dR",&dR);
  selectedTree->SetBranchAddress("relISO_a_corr",&relISO_a_corr);
  selectedTree->SetBranchAddress("ele_pt",&ele_pt);
  selectedTree->SetBranchAddress("ele_sameVertex",&ele_sameVertex);
  selectedTree->SetBranchAddress("ele_ip3D_match",&ele_ip3D_match);
  selectedTree->SetBranchAddress("ele_dxy",&ele_dxy);
  selectedTree->SetBranchAddress("ele_dxyError",&ele_dxyError);
  selectedTree->SetBranchAddress("ele_dz",&ele_dz);
  selectedTree->SetBranchAddress("ele_dzError",&ele_dzError);
  selectedTree->SetBranchAddress("ele_phi",&ele_phi);
  selectedTree->SetBranchAddress("ele_pfPhotonIso",&ele_pfPhotonIso);
  selectedTree->SetBranchAddress("ele_pfChargedHadIso",&ele_pfChargedHadIso);
  selectedTree->SetBranchAddress("ele_pfNeutralHadIso",&ele_pfNeutralHadIso);
  selectedTree->SetBranchAddress("mvaEleID_Fall17_noIso_V2_wp90",&mvaEleID_Fall17_noIso_V2_wp90);
  selectedTree->SetBranchAddress("mvaEleID_Fall17_noIso_V2_wp80",&mvaEleID_Fall17_noIso_V2_wp80);
  selectedTree->SetBranchAddress("mvaEleID_Fall17_iso_V2_wp90",&mvaEleID_Fall17_iso_V2_wp90);
  selectedTree->SetBranchAddress("mvaEleID_Fall17_iso_V2_wp80",&mvaEleID_Fall17_iso_V2_wp80);
  selectedTree->SetBranchAddress("mvaEleID_Fall17_noIso_V2_wpLoose_unsopported",&mvaEleID_Fall17_noIso_V2_wpLoose_unsopported);
  selectedTree->SetBranchAddress("mvaEleID_Fall17_iso_V2_wpHZZ_unsopported",&mvaEleID_Fall17_iso_V2_wpHZZ_unsopported);
  selectedTree->SetBranchAddress("ele_eta",&ele_eta);
  selectedTree->SetBranchAddress("ElectronMVAEstimatorRun2Fall17NoIsoV2Values",&ElectronMVAEstimatorRun2Fall17NoIsoV2Values);
  selectedTree->SetBranchAddress("ElectronMVAEstimatorRun2Fall17IsoV2Values",&ElectronMVAEstimatorRun2Fall17IsoV2Values);
  selectedTree->SetBranchAddress("ele_mother",&ele_mother);
  TBranch *BranchTagPt=selectedTree->Branch("TagPt",&TagPt);
  TBranch *BranchProbePt=selectedTree->Branch("ProbePt",&ProbePt);
  TBranch *BranchTagDxy=selectedTree->Branch("TagDxy",&TagDxy);
  TBranch *BranchProbeDxy=selectedTree->Branch("ProbeDxy",&ProbeDxy);
  TBranch *BranchTagDxyError=selectedTree->Branch("TagDxyError",&TagDxyError);
  TBranch *BranchProbeDxyError=selectedTree->Branch("ProbeDxyError",&ProbeDxyError);
  TBranch *BranchTagDz=selectedTree->Branch("TagDz",&TagDz);
  TBranch *BranchProbeDz=selectedTree->Branch("ProbeDz",&ProbeDz);
  TBranch *BranchTagDzError=selectedTree->Branch("TagDzError",&TagDzError);
  TBranch *BranchProbeDzError=selectedTree->Branch("ProbeDzError",&ProbeDzError);
  TBranch *BranchTagEta=selectedTree->Branch("TagEta",&TagEta);
  TBranch *BranchProbeEta=selectedTree->Branch("ProbeEta",&ProbeEta);
  TBranch *BranchTagPhi=selectedTree->Branch("TagPhi",&TagPhi);
  TBranch *BranchProbePhi=selectedTree->Branch("ProbePhi",&ProbePhi);
  TBranch *BranchTagRelISO=selectedTree->Branch("TagRelISO",&TagRelISO);
  TBranch *BranchProbeRelISO=selectedTree->Branch("ProbeRelISO",&ProbeRelISO);
  TBranch *BranchProbeChIso=selectedTree->Branch("ProbeChIso",&ProbeChIso);
  TBranch *BranchProbeNeuIso=selectedTree->Branch("ProbeNeuIso",&ProbeNeuIso);
  TBranch *BranchProbePhoIso=selectedTree->Branch("ProbePhoIso",&ProbePhoIso);
  TBranch *BranchProbePass=selectedTree->Branch("ProbePass",&ProbePass);
  TBranch *BranchProbeIndex=selectedTree->Branch("ProbeIndex",&ProbeIndex);
  TBranch *BranchTagIndex=selectedTree->Branch("TagIndex",&TagIndex);
  TBranch *BranchTagMVA=selectedTree->Branch("TagMVA",&TagMVA);
  
  Long64_t nentries = selectedTree->GetEntries();
  for(Long64_t i=0;i<nentries;i++) //Tag-Probe pairs, Probe-Tag pairs and Tag-Tag pairs in which first is selected as tag
  {
    selectedTree->GetEntry(i);
    // Define tag and probe tag cuts and pair selection (Check all places where this cut is written. Change needed in 4 lines) 2/4
    if((ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(0)>0.95 && ele_pt->at(0)>7 && (fabs(ele_eta->at(0))<1.44 || fabs(ele_eta->at(0))>1.57)) && (ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(1)>0.95 && ele_pt->at(1)>7 &&(fabs(ele_eta->at(1))<1.44 || fabs(ele_eta->at(1))>1.57)))
    TagIndex=0;
    else
    for(j=0;j<2;j++)
    {
      if(ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(j)>0.95 && ele_pt->at(j)>7 && (fabs(ele_eta->at(j))<1.44 || fabs(ele_eta->at(j))>1.57)) // Define tag and probe tag cuts and pair selection (Check all places where this cut is written. Change needed in 4 lines) 3/4
      TagIndex=j;
    }
    ProbeIndex=abs(TagIndex-1);

    TagPt=ele_pt->at(TagIndex);
    ProbePt=ele_pt->at(ProbeIndex);
    TagDxy=ele_dxy->at(TagIndex);
    ProbeDxy=ele_dxy->at(ProbeIndex);
    TagDxyError=ele_dxyError->at(TagIndex);
    ProbeDxyError=ele_dxyError->at(ProbeIndex);
    TagDz=ele_dz->at(TagIndex);
    ProbeDz=ele_dz->at(ProbeIndex);
    TagDzError=ele_dzError->at(TagIndex);
    ProbeDzError=ele_dzError->at(ProbeIndex);
    TagEta=ele_eta->at(TagIndex);
    ProbeEta=ele_eta->at(ProbeIndex);
    TagPhi=ele_phi->at(TagIndex);
    ProbePhi=ele_phi->at(ProbeIndex);
    TagRelISO=relISO_a_corr->at(TagIndex);
    ProbeRelISO=relISO_a_corr->at(ProbeIndex);
    ProbeChIso=ele_pfChargedHadIso->at(ProbeIndex);
    ProbeNeuIso=ele_pfNeutralHadIso->at(ProbeIndex);
    ProbePhoIso=ele_pfPhotonIso->at(ProbeIndex);
    TagMVA=ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(TagIndex);
    

    // Define probe ID cut/selection (Check all places where this cut is written. Change needed in 2 lines) 1/2

    //if(mvaEleID_Fall17_noIso_V2_wp80->at(ProbeIndex)==1)
    //if(mvaEleID_Fall17_iso_V2_wp80->at(ProbeIndex)==1)
    //if(mvaEleID_Fall17_iso_V2_wpHZZ_unsopported->at(ProbeIndex)==1)
    if(mvaEleID_Fall17_noIso_V2_wpLoose_unsopported->at(ProbeIndex)==1) //ID to measure efficiency

    ProbePass=1;
    else
    ProbePass=0;
    BranchTagPt->Fill();
    BranchTagEta->Fill();
    BranchTagDxy->Fill();
    BranchTagDxyError->Fill();
    BranchTagDz->Fill();
    BranchTagDzError->Fill();
    BranchTagEta->Fill();
    BranchTagPhi->Fill();
    BranchTagRelISO->Fill();
    BranchProbeEta->Fill();
    BranchProbePt->Fill();
    BranchProbeDxy->Fill();
    BranchProbeDxyError->Fill();
    BranchProbeDz->Fill();
    BranchProbeDzError->Fill();
    BranchProbePhi->Fill();
    BranchProbeRelISO->Fill();
    BranchProbeChIso->Fill();
    BranchProbeNeuIso->Fill();
    BranchProbePhoIso->Fill();
    BranchProbePass->Fill();
    BranchTagIndex->Fill();
    BranchProbeIndex->Fill();
    BranchTagMVA->Fill();
  }
  // Define tag and probe tag cuts and pair selection (Check all places where this cut is written. Change needed in 4 lines) 4/4
  TTree* TempTree = originalTree->CopyTree("ele_mother[0]==443 && ele_mother[1]==443 && (ElectronMVAEstimatorRun2Fall17NoIsoV2Values[0]>0.95 && ele_pt[0]>7 && (fabs(ele_eta[0])<1.44 || fabs(ele_eta[0])>1.57)) && (ElectronMVAEstimatorRun2Fall17NoIsoV2Values[1]>0.95 && ele_pt[1]>7 && (fabs(ele_eta[1])<1.44 || fabs(ele_eta[1])>1.57))");
  TempTree->SetBranchAddress("Diele_mass",&Diele_mass);
  TempTree->SetBranchAddress("Diele_pt",&Diele_pt);
  TempTree->SetBranchAddress("dR",&dR);
  TempTree->SetBranchAddress("relISO_a_corr",&relISO_a_corr);
  TempTree->SetBranchAddress("ele_pt",&ele_pt);
  TempTree->SetBranchAddress("ele_sameVertex",&ele_sameVertex);
  TempTree->SetBranchAddress("ele_ip3D_match",&ele_ip3D_match);
  TempTree->SetBranchAddress("ele_dxy",&ele_dxy);
  TempTree->SetBranchAddress("ele_dxyError",&ele_dxyError);
  TempTree->SetBranchAddress("ele_dz",&ele_dz);
  TempTree->SetBranchAddress("ele_dzError",&ele_dzError);
  TempTree->SetBranchAddress("ele_phi",&ele_phi);
  TempTree->SetBranchAddress("ele_pfPhotonIso",&ele_pfPhotonIso);
  TempTree->SetBranchAddress("ele_pfChargedHadIso",&ele_pfChargedHadIso);
  TempTree->SetBranchAddress("ele_pfNeutralHadIso",&ele_pfNeutralHadIso);
  TempTree->SetBranchAddress("mvaEleID_Fall17_noIso_V2_wp90",&mvaEleID_Fall17_noIso_V2_wp90);
  TempTree->SetBranchAddress("mvaEleID_Fall17_noIso_V2_wp80",&mvaEleID_Fall17_noIso_V2_wp80);
  TempTree->SetBranchAddress("mvaEleID_Fall17_iso_V2_wp90",&mvaEleID_Fall17_iso_V2_wp90);
  TempTree->SetBranchAddress("mvaEleID_Fall17_iso_V2_wp80",&mvaEleID_Fall17_iso_V2_wp80);
  TempTree->SetBranchAddress("mvaEleID_Fall17_noIso_V2_wpLoose_unsopported",&mvaEleID_Fall17_noIso_V2_wpLoose_unsopported);
  TempTree->SetBranchAddress("mvaEleID_Fall17_iso_V2_wpHZZ_unsopported",&mvaEleID_Fall17_iso_V2_wpHZZ_unsopported);
  TempTree->SetBranchAddress("ele_eta",&ele_eta);
  TempTree->SetBranchAddress("ElectronMVAEstimatorRun2Fall17NoIsoV2Values",&ElectronMVAEstimatorRun2Fall17NoIsoV2Values);
  TempTree->SetBranchAddress("ElectronMVAEstimatorRun2Fall17IsoV2Values",&ElectronMVAEstimatorRun2Fall17IsoV2Values);
  TempTree->SetBranchAddress("ele_mother",&ele_mother);
  TBranch *BranchTagPt_=TempTree->Branch("TagPt",&TagPt);
  TBranch *BranchProbePt_=TempTree->Branch("ProbePt",&ProbePt);
  TBranch *BranchTagDxy_=TempTree->Branch("TagDxy",&TagDxy);
  TBranch *BranchProbeDxy_=TempTree->Branch("ProbeDxy",&ProbeDxy);
  TBranch *BranchTagDxyError_=TempTree->Branch("TagDxyError",&TagDxyError);
  TBranch *BranchProbeDxyError_=TempTree->Branch("ProbeDxyError",&ProbeDxyError);
  TBranch *BranchTagDz_=TempTree->Branch("TagDz",&TagDz);
  TBranch *BranchProbeDz_=TempTree->Branch("ProbeDz",&ProbeDz);
  TBranch *BranchTagDzError_=TempTree->Branch("TagDzError",&TagDzError);
  TBranch *BranchProbeDzError_=TempTree->Branch("ProbeDzError",&ProbeDzError);
  TBranch *BranchTagEta_=TempTree->Branch("TagEta",&TagEta);
  TBranch *BranchProbeEta_=TempTree->Branch("ProbeEta",&ProbeEta);
  TBranch *BranchTagPhi_=TempTree->Branch("TagPhi",&TagPhi);
  TBranch *BranchProbePhi_=TempTree->Branch("ProbePhi",&ProbePhi);
  TBranch *BranchTagRelISO_=TempTree->Branch("TagRelISO",&TagRelISO);
  TBranch *BranchProbeRelISO_=TempTree->Branch("ProbeRelISO",&ProbeRelISO);
  TBranch *BranchProbeChIso_=TempTree->Branch("ProbeChIso",&ProbeChIso);
  TBranch *BranchProbeNeuIso_=TempTree->Branch("ProbeNeuIso",&ProbeNeuIso);
  TBranch *BranchProbePhoIso_=TempTree->Branch("ProbePhoIso",&ProbePhoIso);
  TBranch *BranchProbePass_=TempTree->Branch("ProbePass",&ProbePass);
  TBranch *BranchProbeIndex_=TempTree->Branch("ProbeIndex",&ProbeIndex);
  TBranch *BranchTagIndex_=TempTree->Branch("TagIndex",&TagIndex);
  TBranch *BranchTagMVA_=TempTree->Branch("TagMVA",&TagMVA);
  nentries = TempTree->GetEntries();
  for(Long64_t i=0;i<nentries;i++) //Tag-Tag pairs in which second is selected as tag
  {
    TempTree->GetEntry(i);
    TagIndex=1;
    ProbeIndex=0;
    TagPt=ele_pt->at(TagIndex);
    ProbePt=ele_pt->at(ProbeIndex);
    TagDxy=ele_dxy->at(TagIndex);
    ProbeDxy=ele_dxy->at(ProbeIndex);
    TagDxyError=ele_dxyError->at(TagIndex);
    ProbeDxyError=ele_dxyError->at(ProbeIndex);
    TagDz=ele_dz->at(TagIndex);
    ProbeDz=ele_dz->at(ProbeIndex);
    TagDzError=ele_dzError->at(TagIndex);
    ProbeDzError=ele_dzError->at(ProbeIndex);
    TagEta=ele_eta->at(TagIndex);
    ProbeEta=ele_eta->at(ProbeIndex);
    TagPhi=ele_phi->at(TagIndex);
    ProbePhi=ele_phi->at(ProbeIndex);
    TagRelISO=relISO_a_corr->at(TagIndex);
    ProbeRelISO=relISO_a_corr->at(ProbeIndex);
    ProbeChIso=ele_pfChargedHadIso->at(ProbeIndex);
    ProbeNeuIso=ele_pfNeutralHadIso->at(ProbeIndex);
    ProbePhoIso=ele_pfPhotonIso->at(ProbeIndex);
    TagMVA=ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(TagIndex);
    
    // Define probe ID cut/selection (Check all places where this cut is written. Change needed in 2 lines) 2/2

    //if(mvaEleID_Fall17_noIso_V2_wp80->at(ProbeIndex)==1)
    //if(mvaEleID_Fall17_iso_V2_wp80->at(ProbeIndex)==1)
    //if(mvaEleID_Fall17_iso_V2_wpHZZ_unsopported->at(ProbeIndex)==1)
    if(mvaEleID_Fall17_noIso_V2_wpLoose_unsopported->at(ProbeIndex)==1) //ID to measure efficiency

    ProbePass=1;
    else
    ProbePass=0;
    BranchTagPt_->Fill();
    BranchTagDxy_->Fill();
    BranchTagDxyError_->Fill();
    BranchTagDz_->Fill();
    BranchTagDzError_->Fill();
    BranchTagEta_->Fill();
    BranchTagPhi_->Fill();
    BranchTagRelISO_->Fill();
    BranchProbeEta_->Fill();
    BranchProbePhi_->Fill();
    BranchProbePt_->Fill();
    BranchProbeDxy_->Fill();
    BranchProbeDxyError_->Fill();
    BranchProbeDz_->Fill();
    BranchProbeDzError_->Fill();
    BranchProbeRelISO_->Fill();
    BranchProbePass_->Fill();
    BranchTagIndex_->Fill();
    BranchProbeIndex_->Fill();
    BranchTagMVA_->Fill();
  }

  TList *list = new TList;
  list->Add(selectedTree);
  list->Add(TempTree);
  TTree *FinalTree = TTree::MergeTrees(list);
  FinalTree->Write();
}
