#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"

using namespace RooFit;

std::pair<float, float> MakeFitDataAltBgk(TTree* tree, std::tuple<float,float,float,float,float,float,float,float,float> dInfo, bool DoPlot, std::string FigName, bool fixedMass)
{
  RooRealVar Diele_mass("Diele_mass", "Dielectron mass", 2,4);
  RooArgSet ntupleVarSet(Diele_mass);
  RooDataSet DataSet("data", "data set", tree, ntupleVarSet);

  RooRealVar dcbMean("dcbMean","dcbMean",std::get<0>(dInfo),std::get<0>(dInfo)*0.995,std::get<0>(dInfo)*1.015);//5

  if(fixedMass){
    dcbMean.setVal(std::get<0>(dInfo)); // Fix JPsi mass
    dcbMean.setConstant(kTRUE); //Fix Jpsi mass
  }


  RooRealVar dcbSigma("dcbSigma","dcbSigma",std::get<1>(dInfo),std::get<1>(dInfo)*0.985,std::get<1>(dInfo)*1.025);
  RooRealVar dcbAlphaL("dcbAlphaL","dcbAlphaL",std::get<2>(dInfo),std::get<2>(dInfo)*0.95,std::get<2>(dInfo)*1.05);
  RooRealVar dcbNL("dcbNL","dcbNL",std::get<3>(dInfo),std::get<3>(dInfo)*0.95,std::get<3>(dInfo)*1.05);
  RooRealVar dcbAlphaR("dcbAlphaR","dcbAlphaR",std::get<4>(dInfo),std::get<4>(dInfo)*0.95,std::get<4>(dInfo)*1.05);
  RooRealVar dcbNR("dcbNR","dcbNR",std::get<5>(dInfo),std::get<5>(dInfo)*0.95,std::get<5>(dInfo)*1.05);
  RooCrystalBall dcb("dcb","dcb",Diele_mass,dcbMean,dcbSigma,dcbAlphaL,dcbNL,dcbAlphaR,dcbNR);
  RooRealVar sig_yield("sig_yield", "yield of signal peak", 770, 300, 1000000);



  RooRealVar sl("sl", "sl", -0.5, -10, 10);
  RooExponential bkg_poly("background","background", Diele_mass, sl);

  RooRealVar bkg_yield("bkg_yield", "yield of background", 50000, 0, 10000000);

  RooArgList shapes, yields;
  shapes.add(bkg_poly);
  yields.add(bkg_yield);
  shapes.add(dcb);
  yields.add(sig_yield);
  RooAddPdf totalPdf("totalPdf", "sum of signal and background PDF's", shapes, yields);

  totalPdf.fitTo(DataSet, Extended());

  if(DoPlot)
  {
  TCanvas *canvas = new TCanvas("cs_data","cs_data",1600,800);
  RooPlot* massframe = Diele_mass.frame(Title("alt. bgk +1"));
  DataSet.plotOn(massframe);
  totalPdf.plotOn(massframe,Components(bkg_poly), LineColor(kRed));
  totalPdf.plotOn(massframe,Components(dcb), LineColor(kMagenta));
  totalPdf.plotOn(massframe, LineColor(kBlue));
  totalPdf.paramOn(massframe, Layout(0.75));
  massframe->Draw();
  const char *c = FigName.c_str();
  canvas->SaveAs(c);
  }

  return std::make_pair(sig_yield.getVal(), sig_yield.getError());
}
