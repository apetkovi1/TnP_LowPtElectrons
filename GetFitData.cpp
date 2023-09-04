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
#include <tuple>

using namespace RooFit;

std::tuple<float,float,float,float,float,float,float,float,float> GetFitData(TTree* tree, std::pair<float,float> mcInfo)
{
  RooRealVar Diele_mass("Diele_mass", "Dielectron mass", 2,4);
  RooArgSet ntupleVarSet(Diele_mass);
  RooDataSet DataSet("data", "data set", tree, ntupleVarSet);

  RooRealVar dcbMean("dcbMean","dcbMean",std::get<0>(mcInfo),std::get<0>(mcInfo)*0.995,std::get<0>(mcInfo)*1.005);
  RooRealVar dcbSigma("dcbSigma","dcbSigma",std::get<1>(mcInfo),std::get<1>(mcInfo)*0.95,std::get<1>(mcInfo)*1.1);
  dcbMean.setVal(std::get<0>(mcInfo)); // Fix JPsi mass
  dcbMean.setConstant(kTRUE); //Fix Jpsi mass
  RooRealVar dcbAlphaL("dcbAlphaL","dcbAlphaL",0.57, 0.3, 10.);
  RooRealVar dcbNL("dcbNL","dcbNL",30., 5., 60.);
  RooRealVar dcbAlphaR("dcbAlphaR","dcbAlphaR",1.1, 0.3, 10.);
  RooRealVar dcbNR("dcbNR","dcbNR",30., 5., 60.);
  //  dcbMean.setVal(mcMean); //3.096916); // Fix JPsi mass
  //  dcbMean.setConstant(kTRUE); //Fix Jpsi mass
  RooCrystalBall dcb("dcb","dcb",Diele_mass,dcbMean,dcbSigma,dcbAlphaL,dcbNL,dcbAlphaR,dcbNR);
  RooRealVar sig_yield("sig_yield", "yield of signal peak", 770, 0, 1000000);

  RooRealVar a0("a0", "a0", -0.5, -0.75, 1);
  RooRealVar a1("a1", "a1", 0.08, 0.05, 0.15);
  RooRealVar a2("a2", "a2", 0.01, 0.003, 1);
  RooChebychev bkg_poly("background","background", Diele_mass, RooArgList(a0,a1,a2));
  RooRealVar bkg_yield("bkg_yield", "yield of background", 50000, 0, 10000000);

  RooArgList shapes, yields;
  shapes.add(bkg_poly);
  yields.add(bkg_yield);
  shapes.add(dcb);
  yields.add(sig_yield);
  RooAddPdf totalPdf("totalPdf", "sum of signal and background PDF's", shapes, yields);

  totalPdf.fitTo(DataSet, Extended());


  return std::make_tuple(dcbMean.getVal(),dcbSigma.getVal(),dcbAlphaL.getVal(),dcbNL.getVal(),dcbAlphaR.getVal(),dcbNR.getVal(),a0.getVal(),a1.getVal(),a2.getVal());
}
