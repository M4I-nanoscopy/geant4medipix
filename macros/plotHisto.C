{
  gROOT->Reset();
  gStyle->SetHistFillColor(kAzure -9 ); 
  gStyle->SetHistLineColor(kAzure -7 ); 
  
  // Draw histos filled by Geant4 simulation 
  //   
  //TFile f = TFile("gamma.root");
  TFile f = TFile("electron.root");
  //TFile f = TFile("proton.root");  
  //TFile f = TFile("alpha.root");  

  TCanvas* c1 = new TCanvas("c1", "  ");
  TH1D* hist1 = (TH1D*)f.Get("1;1");
  c1->cd();
  hist1->Draw("F HIST"); 
  hist1->UseCurrentStyle();
  hist1->GetXaxis()->SetTitle("Depth(um)");
  hist1->GetYaxis()->SetTitle("Counts");
  c1->Update();
  
  TCanvas* c2 = new TCanvas("c2", "  ");
  TH1D* hist2 = (TH1D*)f.Get("2;1");
  hist2->Draw("HIST");
  hist1->GetXaxis()->SetTitle("Energy(keV)");
  hist1->GetYaxis()->SetTitle("Counts");
  hist2->UseCurrentStyle();
  c2->Update();
  
  TCanvas* c3 = new TCanvas("c3", "  ");
  TH1D* hist3 = (TH1D*)f.Get("3;1");
  hist3->Draw("HIST");
  hist1->GetXaxis()->SetTitle("Energy(keV)");
  hist1->GetYaxis()->SetTitle("Counts");
  hist3->UseCurrentStyle();
  c3->Update();
}  
