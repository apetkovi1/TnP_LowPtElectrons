
TEfficiency* Efficiency(std::vector<std::pair<float,float>> YieldErrorPass, std::vector<std::pair<float,float>> YieldErrorFail, float bins[], int NumberOfBins)
{

  TH1F* hist_pass = new TH1F("hist pass", "hist pass",NumberOfBins, bins);
  TH1F* hist_fail = new TH1F("hist fail", "hist fail",NumberOfBins, bins);
  TH1F* hist_total = new TH1F("hist total", "hist total",NumberOfBins, bins);
  for (int i = 0; i < NumberOfBins; i++)
   {
       hist_pass->SetBinContent(i+1, YieldErrorPass.at(i).first);
       hist_pass->SetBinError(i+1,YieldErrorPass.at(i).second);
       hist_fail->SetBinContent(i+1,YieldErrorFail.at(i).first);
       hist_fail->SetBinError(i+1,YieldErrorFail.at(i).second);
   }

   hist_total->Add(hist_pass);
   hist_total->Add(hist_fail);

   TEfficiency* Eff= new TEfficiency(*hist_pass,*hist_total);

  /*Eff_data->SetLineColor(kBlue);
  Eff_MC->SetLineColor(kRed);
  auto legend = new TLegend(0.9,0.9,1.0,1.0);
  auto mg  = new TMultiGraph();
  TCanvas* oi = new TCanvas();
  oi->cd();
  Eff_data->Draw();
  gPad->Update();
  auto graph_data = Eff_data->GetPaintedGraph();
  mg->Add(graph_data);
  Eff_MC->Draw();
  gPad->Update();
  auto graph_MC = Eff_MC->GetPaintedGraph();
  mg->Add(graph_MC);
  legend->AddEntry(graph_data,"data");
  legend->AddEntry(graph_MC,"MC");
  mg->Draw("ap");
  mg->GetXaxis()->SetTitle("pT/GeV");
  mg->GetYaxis()->SetTitle("Efficiency");
  legend->Draw();
  oi->SaveAs("plots/Efficiency.pdf");*/

   return Eff;
}
