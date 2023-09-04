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

std::pair<float, float> GetFitMC(TTree* tree)
{
  RooRealVar Diele_mass("Diele_mass", "Dielectron mass", 2,4);
  RooArgSet ntupleVarSet(Diele_mass);
  RooDataSet DataSet("data", "data set", tree, ntupleVarSet);

  RooRealVar dcbMean("dcbMean","dcbMean",3.1, 3.0,3.15);
  RooRealVar dcbSigma("dcbSigma","dcbSigma",0.085, 0.053, 1.);
  RooRealVar dcbAlphaL("dcbAlphaL","dcbAlphaL",0.57, 0.1, 10.);
  RooRealVar dcbNL("dcbNL","dcbNL",2.8, 0, 60.);
  RooRealVar dcbAlphaR("dcbAlphaR","dcbAlphaR",1.1, 0.1, 10.);
  RooRealVar dcbNR("dcbNR","dcbNR",1.7, 0, 60.);
  RooCrystalBall dcb("dcb","dcb",Diele_mass,dcbMean,dcbSigma,dcbAlphaL,dcbNL,dcbAlphaR,dcbNR);
  RooRealVar sig_yield("sig_yield", "yield of signal peak", 1000, 0, 1000000);


  RooArgList shapes, yields;
  shapes.add(dcb);
  yields.add(sig_yield);
  RooAddPdf totalPdf("totalPdf", "sum of signal and background PDF's", shapes, yields);

  totalPdf.fitTo(DataSet, Extended());


  return std::make_pair(dcbMean.getVal(),dcbSigma.getVal());
}
