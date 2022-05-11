
void Efficiency(std::vector<std::pair<float,float>> YieldErrorDATAPass, std::vector<std::pair<float,float>> YieldErrorDataFail,
                std::vector<std::pair<float,float>> YieldErrorMCPass, std::vector<std::pair<float,float>> YieldErrorMCFail, float bins[], int NumberOfBins)
{

  TH1F* hist_pass_data = new TH1F("hist pass data", "hist pass data",NumberOfBins, bins);
  TH1F* hist_fail_data = new TH1F("hist fail data", "hist fail data",NumberOfBins, bins);
  TH1F* hist_pass_MC = new TH1F("hist pass MC", "hist pass MC",NumberOfBins, bins);
  TH1F* hist_fail_MC = new TH1F("hist fail MC", "hist fail MC",NumberOfBins, bins);
  TH1F* hist_total_data = new TH1F("hist total data", "hist total data",NumberOfBins, bins);
  TH1F* hist_total_MC = new TH1F("hist total MC", "hist total MC",NumberOfBins, bins);
  for (int i = 0; i < NumberOfBins; i++)
   {
       hist_pass_data->SetBinContent(i+1, YieldErrorDATAPass.at(i).first);
       hist_pass_data->SetBinError(i+1,YieldErrorDATAPass.at(i).second);
       hist_fail_data->SetBinContent(i+1,YieldErrorDataFail.at(i).first);
       hist_fail_data->SetBinError(i+1,YieldErrorDataFail.at(i).second);

       hist_pass_MC->SetBinContent(i+1,YieldErrorMCPass.at(i).first);
       hist_pass_MC->SetBinError(i+1,YieldErrorMCPass.at(i).second);
       hist_fail_MC->SetBinContent(i+1,YieldErrorMCFail.at(i).first);
       hist_fail_MC->SetBinError(i+1,YieldErrorMCFail.at(i).second);
   }

   hist_total_data->Add(hist_pass_data);
   hist_total_data->Add(hist_fail_data);

   hist_total_MC->Add(hist_pass_MC);
   hist_total_MC->Add(hist_fail_MC);

   TEfficiency* Eff_data= new TEfficiency(*hist_pass_data,*hist_total_data);
   TEfficiency* Eff_MC= new TEfficiency(*hist_pass_MC,*hist_total_MC);

  Eff_data->SetLineColor(kBlue);
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
  oi->SaveAs("plots/Efficiency.pdf");
}
