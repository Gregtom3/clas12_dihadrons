
void perform_RooUnfold(int usable_bins, std::vector<int> true_bins, std::vector<int> reco_bins, std::vector<int> data_bins, const char* method="bayes") {
  double bw = 0.5;
  TFile *f = new TFile("/work/clas12/users/gmat/tmp/unfolded_hists.root", "RECREATE");
  RooUnfoldResponse response(usable_bins, 1.0, usable_bins + 1.0);
  TH1D *true_hist = new TH1D("true_hist", "", usable_bins, 1.0, usable_bins + 1.0);
  TH1D *reco_hist = new TH1D("reco_hist", "", usable_bins, 1.0, usable_bins + 1.0);
  TH1D *data_hist = new TH1D("data_hist", "", usable_bins, 1.0, usable_bins + 1.0);

  cout << true_bins.size() << endl;
  for (int i = 0; i < true_bins.size(); i++) {
    int tb = true_bins[i];
    int rb = reco_bins[i];

    if (tb == 0 || rb == 0) continue;
    if (tb == -1 && rb == -1) continue;

    if (tb != -1) true_hist->Fill(tb + bw);
    if (rb == -1 && tb != -1) response.Miss(tb + bw);
    else if (rb != -1 && tb == -1) {
      response.Fake(rb + bw);
      reco_hist->Fill(rb + bw);
    } else {
      response.Fill(rb + bw, tb + bw);
      reco_hist->Fill(rb + bw);
    }
  }
    cout << true_hist->GetMean() << endl;
  for (int db : data_bins) {
    if (db == 0) continue;
    data_hist->Fill(db + bw);
  }

  TH1D *hUnfold;

  if (strcmp(method, "bayes") == 0) {
    RooUnfoldBayes *uf = new RooUnfoldBayes(&response, data_hist, 25);
    uf->HandleFakes(1);
    hUnfold = (TH1D*) uf->Hunfold();
  } else if (strcmp(method, "svd") == 0) {
    RooUnfoldSvd *uf = new RooUnfoldSvd(&response, data_hist, 20);
    hUnfold = (TH1D*) uf->Hunfold();
  }


  reco_hist->Write();
  true_hist->Write();
  data_hist->Write();
  hUnfold->Write();
  f->Close();
}