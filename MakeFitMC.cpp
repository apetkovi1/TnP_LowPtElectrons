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

std::pair<float, float> MakeFitMC(TTree* tree, bool DoPlot, std::string FigName)
{
  RooRealVar Diele_mass("Diele_mass", "Dielectron mass", 2,4);
  RooArgSet ntupleVarSet(Diele_mass);
  RooDataSet DataSet("data", "data set", tree, ntupleVarSet);

  RooRealVar* dcbMean  = new RooRealVar("dcbMean","dcbMean",3.1, 3.0,3.15);
  RooRealVar* dcbSigma = new RooRealVar("dcbSigma","dcbSigma",0.05, 0, 1.);
  RooRealVar* dcbAlphaL = new RooRealVar("dcbAlphaL","dcbAlphaL",0.57, 0, 10.);
  RooRealVar* dcbNL = new RooRealVar("dcbNL","dcbNL",2.8, 0, 10.);
  RooRealVar* dcbAlphaR = new RooRealVar("dcbAlphaR","dcbAlphaR",1.1, 0, 10.);
  RooRealVar* dcbNR = new RooRealVar("dcbNR","dcbNR",1.7, 0, 10.);
  RooCrystalBall* dcb = new RooCrystalBall("dcb","dcb",Diele_mass,*dcbMean,*dcbSigma, *dcbAlphaL, *dcbNL, *dcbAlphaR, *dcbNR);
  RooRealVar sig_yield("sig_yield", "yield of signal peak", 1000, 0, 1000000);


  RooArgList shapes, yields;
  shapes.add(*dcb);
  yields.add(sig_yield);
  RooAddPdf totalPdf("totalPdf", "sum of signal and background PDF's", shapes, yields);

  totalPdf.fitTo(DataSet, Extended());

  if(DoPlot)
  {
  TCanvas *canvas = new TCanvas("cs_data","cs_data",1600,800);
  RooPlot* massframe = Diele_mass.frame();
  DataSet.plotOn(massframe);
  totalPdf.plotOn(massframe, LineColor(kBlue));
  totalPdf.paramOn(massframe, Layout(0.75));
  massframe->Draw();
  const char *c = FigName.c_str();
  canvas->SaveAs(c);
  }

  return std::make_pair(sig_yield.getVal(), sig_yield.getError());
}
