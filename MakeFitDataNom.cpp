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

std::pair<float, float> MakeFitDataNom(TTree* tree, std::tuple<float,float,float,float,float,float,float,float,float> dInfo, bool DoPlot, std::string FigName, bool fixedMass)
{
  RooRealVar Diele_mass("Diele_mass", "Dielectron mass", 2,4);
  RooArgSet ntupleVarSet(Diele_mass);
  RooDataSet DataSet("data", "data set", tree, ntupleVarSet);

  RooRealVar dcbMean("dcbMean","dcbMean",std::get<0>(dInfo),std::get<0>(dInfo)*0.995,std::get<0>(dInfo)*1.005);
  // RooRealVar dcbSigma("dcbSigma","dcbSigma",std::get<1>(dInfo),0.03,std::get<1>(dInfo)*1.2);
  RooRealVar dcbSigma("dcbSigma","dcbSigma",std::get<1>(dInfo),std::get<1>(dInfo)*0.955,std::get<1>(dInfo)*1.15);
  //RooRealVar dcbAlphaL("dcbAlphaL","dcbAlphaL",std::get<2>(dInfo),2,std::get<2>(dInfo)*1.2);
  if(fixedMass){
    dcbMean.setVal(std::get<0>(dInfo)); // Fix JPsi mass
    dcbMean.setConstant(kTRUE); //Fix Jpsi mass
  }
  RooRealVar dcbAlphaL("dcbAlphaL","dcbAlphaL",std::get<2>(dInfo),std::get<2>(dInfo)*0.8,std::get<2>(dInfo)*1.2);
  RooRealVar dcbNL("dcbNL","dcbNL",std::get<3>(dInfo),std::get<3>(dInfo)*0.9,std::get<3>(dInfo)*1.2);
  // RooRealVar dcbAlphaR("dcbAlphaR","dcbAlphaR",std::get<4>(dInfo),2,std::get<4>(dInfo)*1.2);
  RooRealVar dcbAlphaR("dcbAlphaR","dcbAlphaR",std::get<4>(dInfo),std::get<4>(dInfo)*0.8,std::get<4>(dInfo)*1.2);
  RooRealVar dcbNR("dcbNR","dcbNR",std::get<5>(dInfo),std::get<5>(dInfo)*0.9,std::get<5>(dInfo)*1.2);
  //  dcbMean->setVal(3.096916); // Fix JPsi mass
  //  dcbMean->setConstant(kTRUE); //Fix Jpsi mass
  RooCrystalBall dcb("dcb","dcb",Diele_mass,dcbMean,dcbSigma,dcbAlphaL,dcbNL,dcbAlphaR,dcbNR);
  RooRealVar sig_yield("sig_yield", "yield of signal peak", 770, 0, 1000000);

  RooRealVar a0("a0", "a0", std::get<6>(dInfo),-0.75,1);
  RooRealVar a1("a1", "a1", std::get<7>(dInfo),0.05,0.15);
  RooRealVar a2("a2", "a2", std::get<8>(dInfo),0.003,1);

  /*

 RooRealVar a0("a0", "a0", std::get<6>(dInfo),std::get<6>(dInfo)*0.9,std::get<6>(dInfo)*1.1);
  RooRealVar a1("a1", "a1", std::get<7>(dInfo),std::get<7>(dInfo)*0.9,std::get<7>(dInfo)*1.1);
RooRealVar a2("a2", "a2", std::get<8>(dInfo),std::get<8>(dInfo)*0.9,std::get<8>(dInfo)*1.1);

  RooRealVar dcbMean("dcbMean","dcbMean",std::get<0>(dInfo));
  RooRealVar dcbSigma("dcbSigma","dcbSigma",std::get<1>(dInfo));
  RooRealVar dcbAlphaL("dcbAlphaL","dcbAlphaL",std::get<2>(dInfo));
  RooRealVar dcbNL("dcbNL","dcbNL",std::get<3>(dInfo));
  RooRealVar dcbAlphaR("dcbAlphaR","dcbAlphaR",std::get<4>(dInfo));
  RooRealVar dcbNR("dcbNR","dcbNR",std::get<5>(dInfo));
  //  dcbMean->setVal(3.096916); // Fix JPsi mass
  //  dcbMean->setConstant(kTRUE); //Fix Jpsi mass
  RooCrystalBall dcb("dcb","dcb",Diele_mass,dcbMean,dcbSigma,dcbAlphaL,dcbNL,dcbAlphaR,dcbNR);
  RooRealVar sig_yield("sig_yield", "yield of signal peak", 770, 0, 1000000);

  RooRealVar a0("a0", "a0", std::get<6>(dInfo));
  RooRealVar a1("a1", "a1", std::get<7>(dInfo));
  RooRealVar a2("a2", "a2", std::get<8>(dInfo));
  
  */

  RooChebychev bkg_poly("background","background", Diele_mass, RooArgList(a0,a1,a2));
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
      RooPlot* massframe = Diele_mass.frame(Title("nominal fit"));
      //      RooPlot* massframe = Diele_mass.frame();
      DataSet.plotOn(massframe);
      totalPdf.plotOn(massframe,Components(bkg_poly), LineColor(kRed));
      totalPdf.plotOn(massframe,Components(dcb), LineColor(kMagenta));
      totalPdf.plotOn(massframe, LineColor(kBlue));
      totalPdf.paramOn(massframe, Layout(0.75));
      massframe->Draw();
      const char *c = FigName.c_str();
      canvas->SaveAs(c);
    }
  std::cout << "\n\n\n\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  std::cout << "yield: " << sig_yield.getVal() << "\tits err: " << sig_yield.getError() << std::endl;
  return std::make_pair(sig_yield.getVal(), sig_yield.getError());
}
