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

std::pair<float, float> MakeFitData(TTree* tree, bool DoPlot, std::string FigName)
{
  RooRealVar Diele_mass("Diele_mass", "Dielectron mass", 2,4);
  RooArgSet ntupleVarSet(Diele_mass);
  RooDataSet DataSet("data", "data set", tree, ntupleVarSet);

  RooRealVar* dcbMean  = new RooRealVar("dcbMean","dcbMean",3.1, 3.05,3.16);
  //RooRealVar* dcbMean  = new RooRealVar("dcbMean","dcbMean",3.096916, 3.096916,3.096916);
  RooRealVar* dcbSigma = new RooRealVar("dcbSigma","dcbSigma",0.16, 0.05, 0.30); // 0.05
  RooRealVar* dcbAlphaL = new RooRealVar("dcbAlphaL","dcbAlphaL",2.05, 0.5, 10.); // 0.57
  RooRealVar* dcbNL = new RooRealVar("dcbNL","dcbNL",30., 10., 95.);// min value 20, uzliku 10 un 5, no 30 un 10
  RooRealVar* dcbAlphaR = new RooRealVar("dcbAlphaR","dcbAlphaR",4.4, 0.5, 9.);
  RooRealVar* dcbNR = new RooRealVar("dcbNR","dcbNR",30., 10., 95.); // min value 20
  RooCrystalBall* dcb = new RooCrystalBall("dcb","dcb",Diele_mass,*dcbMean,*dcbSigma, *dcbAlphaL, *dcbNL, *dcbAlphaR, *dcbNR);
  RooRealVar sig_yield("sig_yield", "yield of signal peak", 400, 30, 1000000); //770

  RooRealVar a0("a0", "a0", -0.5, -2, 2);// bija -10 lidz 10 abos
  RooRealVar a1("a1", "a1", -0.03, -0.4, 0.5);
  RooChebychev bkg_poly("background","background", Diele_mass, RooArgList(a0,a1));
  RooRealVar bkg_yield("bkg_yield", "yield of background", 50000, 0, 10000000);

  RooArgList shapes, yields;
  shapes.add(bkg_poly);
  yields.add(bkg_yield);
  shapes.add(*dcb);
  yields.add(sig_yield);
  RooAddPdf totalPdf("totalPdf", "sum of signal and background PDF's", shapes, yields);

  totalPdf.fitTo(DataSet, Extended());

  if(DoPlot)
  {
  TCanvas *canvas = new TCanvas("cs_data","cs_data",1600,800);
  RooPlot* massframe = Diele_mass.frame();
  DataSet.plotOn(massframe);
  totalPdf.plotOn(massframe,Components(bkg_poly), LineColor(kRed));
  totalPdf.plotOn(massframe,Components(*dcb), LineColor(kMagenta));
  totalPdf.plotOn(massframe, LineColor(kBlue));
  totalPdf.paramOn(massframe, Layout(0.75));
  massframe->Draw();
  const char *c = FigName.c_str();
  canvas->SaveAs(c);
  }

  return std::make_pair(sig_yield.getVal(), sig_yield.getError());
}
