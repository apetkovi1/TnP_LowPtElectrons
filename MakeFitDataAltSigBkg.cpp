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

std::pair<float, float> MakeFitDataAltSigBkg(TTree* tree, std::tuple<float,float,float,float,float,float,float,float,float> dInfo, bool DoPlot, std::string FigName,bool fixedMass)
{
  RooRealVar Diele_mass("Diele_mass", "Dielectron mass", 2,4);
  RooArgSet ntupleVarSet(Diele_mass);
  RooDataSet DataSet("data", "data set", tree, ntupleVarSet);

  RooRealVar sigMean("gMean","gMean",std::get<0>(dInfo),std::get<0>(dInfo)*0.985,std::get<0>(dInfo)*1.025);

  if(fixedMass){
    sigMean.setVal(std::get<0>(dInfo)); // Fix JPsi mass
    sigMean.setConstant(kTRUE); //Fix Jpsi mass
  }


  RooRealVar sigSigma("gSigma","gSigma",std::get<1>(dInfo),std::get<1>(dInfo)*0.9,std::get<1>(dInfo)*1.1);
  RooGaussian gauss("gauss","gaussian PDF",Diele_mass,sigMean,sigSigma);


  RooRealVar sl("sl", "slope", 8.1, 0.5,215);
  RooExponential exp("exp", "exp", Diele_mass, sl);
  
  RooGenericPdf step_func("step_func","step_func","@0 < 0.0 ? @1 : 0.0",RooArgSet(Diele_mass, exp));
  RooFFTConvPdf gxe("gxe", "gauss (X) exponential", Diele_mass, gauss, step_func);
  RooRealVar sig_yield("sig_yield", "yield of signal peak", 770, 0, 1000000);

  RooRealVar slope("slope_bkg","slope_bkg", -0.5 , -10 , 10 );
  RooExponential bkg_poly("background", "background", Diele_mass, slope);
  RooRealVar bkg_yield("bkg_yield", "yield of background", 50000, 0, 10000000);

  RooArgList shapes, yields;
  shapes.add(bkg_poly);
  yields.add(bkg_yield);
  shapes.add(gxe);
  yields.add(sig_yield);
  RooAddPdf totalPdf("totalPdf", "sum of signal and background PDF's", shapes, yields);

  totalPdf.fitTo(DataSet, Extended());

  if(DoPlot)
    {
      TCanvas *canvas = new TCanvas("cs_data","cs_data",1600,800);
      RooPlot* massframe = Diele_mass.frame(Title("exp (x) gauss , +altBkg"));
      DataSet.plotOn(massframe);
      totalPdf.plotOn(massframe,Components(bkg_poly), LineColor(kRed));
      totalPdf.plotOn(massframe,Components(gxe), LineColor(kMagenta));
      totalPdf.plotOn(massframe, LineColor(kBlue));
      totalPdf.paramOn(massframe, Layout(0.75));
      massframe->Draw();
      const char *c = FigName.c_str();
      canvas->SaveAs(c);
    }

  return std::make_pair(sig_yield.getVal(), sig_yield.getError());
}
