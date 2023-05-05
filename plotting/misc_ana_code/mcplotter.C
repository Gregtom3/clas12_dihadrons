std::vector<double> linspace(double start, double stop, int num) {
    std::vector<double> bins;
    double step = (stop - start) / (num - 1);

    for (int i = 0; i < num; ++i) {
        bins.push_back(start + i * step);
    }

    return bins;
}

std::vector<double> logspace_bins(double start, double stop, int num) {
    std::vector<double> bins;
    double log_start = std::log10(start);
    double log_stop = std::log10(stop);
    double step = (log_stop - log_start) / (num - 1);

    for (int i = 0; i < num; ++i) {
        double value = log_start + i * step;
        bins.push_back(std::pow(10, value));
    }

    return bins;
}
std::vector<std::vector<std::string>> generate_input_boundaries(const std::vector<std::pair<double, double>>& x_bins, const std::vector<int>& add_middle_cut) {
    std::string low_Q2_cut = "Q2>0.9253 + 3.8613 * x + -9.4351 * x*x + 16.9292 * x*x*x + 17.2332 * x*x*x*x";
    std::string middle_Q2_cut = "0.6361935324019532 + 5.961630973508846*x + 12.028695097029118*x*x";
    std::string high_Q2_cut = "Q2<x*17";
    std::vector<std::vector<std::string>> input_boundaries;

    for (size_t i = 0; i < x_bins.size(); ++i) {
        double x_min = x_bins[i].first;
        double x_max = x_bins[i].second;
        std::vector<std::string> input_boundary;
        if(add_middle_cut[i]==0){
          input_boundary = {
            "x<=" + std::to_string(x_max),
            "x>=" + std::to_string(x_min),
            high_Q2_cut,
            low_Q2_cut
          };
        }else{
          input_boundary = {
              "x<=" + std::to_string(x_max),
              "x>=" + std::to_string(x_min),
              "",
              high_Q2_cut,
              low_Q2_cut
          };
          if (add_middle_cut[i] == -1) {
              input_boundary[2] = "Q2<" + middle_Q2_cut;
          } else if (add_middle_cut[i] == 1) {
              input_boundary[2] = "Q2>" + middle_Q2_cut;
          }
        }
        input_boundaries.push_back(input_boundary);
    }

    return input_boundaries;
}


std::vector<std::vector<std::string>> generate_z_pt_boundaries(const std::vector<double>& z_edges, const std::vector<double>& pT_edges) {
    std::vector<std::vector<std::string>> boundaries;
    int nzbins = z_edges.size() - 1;
    int npTbins = pT_edges.size() - 1;

    for (int iz = 0; iz < nzbins; ++iz) {
        for (int ip = 0; ip < npTbins; ++ip) {
            boundaries.push_back({ "z>=" + std::to_string(z_edges[iz]),
                                    "z<" + std::to_string(z_edges[iz + 1]),
                                    "pTtot>=" + std::to_string(pT_edges[ip]),
                                    "pTtot<" + std::to_string(pT_edges[ip + 1]) });
        }
    }

    return boundaries;
}


std::string create_filter_string(const std::vector<std::string>& boundary) {
    std::string filter_string;
    for (size_t i = 0; i < boundary.size(); ++i) {
        filter_string += boundary[i];
        if (i < boundary.size() - 1) {
            filter_string += " && ";
        }
    }
    return filter_string;
}




std::tuple<TH2D, std::vector<TH2D>, std::vector<TH2D>, std::vector<std::vector<TH2D>>> get_hists(TString infile, std::vector<std::vector<std::string>> xQ2_boundaries, std::vector<std::vector<std::string>> zpT_boundaries, TString var1, TString var2, std::vector<double> var1_bins, std::vector<double> var2_bins, TString var3, std::vector<double> var3_bins, TString var4, std::vector<double> var4_bins) {

  // define cut conditions
  string cut_conditions = "xF1>0&&xF2>0&&z<0.95&&Mx>1.5";
    
  // create RDataFrame and select events that satisfy the cut conditions
  ROOT::RDataFrame df("dihadron", infile.Data());
  auto df_cut = df.Filter(cut_conditions);

  // create main histogram
  auto main_hist = df_cut.Histo2D({TString(Form("h2d_main_%s_%s", var1.Data(), var2.Data())).Data(), TString(Form("%s-%s distribution", var1.Data(), var2.Data())).Data(), (int)(var1_bins.size() - 1), var1_bins.data(), (int)(var2_bins.size() - 1), var2_bins.data()}, string(var1), string(var2));
    
  std::vector<ROOT::RDF::RResultPtr<TH2D>> var1_var2_hists;
  std::vector<ROOT::RDF::RResultPtr<TH2D>> var3_var4_hists;
  std::vector<std::vector<ROOT::RDF::RResultPtr<TH2D>>> var3_var4_subhists;
    
  std::vector<ROOT::RDF::RNode> tmp_df_cuts;
  for (size_t idx = 0; idx < xQ2_boundaries.size(); ++idx) {
    std::cout << (idx + 1) << " of " << xQ2_boundaries.size() << std::endl;
    //auto tmp_df_cut = df_cut.Filter(xQ2_boundaries[idx][0] + "&&" + xQ2_boundaries[idx][1]);
    //tmp_df_cuts.push_back(tmp_df_cut);
        
    // create histogram for var1 and var2
    std::string filter_string = create_filter_string(xQ2_boundaries[idx]);
    auto df_cut_cut = df_cut.Filter(filter_string);
    auto h2d_var1_var2 = df_cut_cut.Histo2D({TString(Form("h2d_%s_%s_%zu", var1.Data(), var2.Data(), idx)).Data(), TString(Form("%s-%s distribution", var1.Data(), var2.Data())).Data(), (int)(var1_bins.size() - 1), var1_bins.data(), (int)(var2_bins.size() - 1), var2_bins.data()}, string(var1), string(var2));
        
    var1_var2_hists.push_back(h2d_var1_var2);

    auto h2d_var3_var4 = df_cut_cut.Histo2D({TString(Form("h2d_%s_%s_%zu", var3.Data(), var4.Data(), idx)).Data(), TString(Form("%s-%s distribution", var3.Data(), var4.Data())).Data(), (int)(var3_bins.size() - 1), var3_bins.data(), (int)(var4_bins.size() - 1), var4_bins.data()}, string(var3), string(var4));
      
    var3_var4_hists.push_back(h2d_var3_var4);
    
    std::vector<ROOT::RDF::RResultPtr<TH2D>> var3_var4_subhists_arr;
    
    for (size_t idxidx = 0; idxidx < zpT_boundaries.size(); ++idxidx) {
      //auto tmptmp_df_cut = tmp_df_cut.Filter(zpT_boundaries[idxidx][0] + "&&" + zpT_boundaries[idxidx][1]);
      //tmp_df_cuts.push_back(tmptmp_df_cut);
      std::string filter_string = create_filter_string(zpT_boundaries[idxidx]);
      auto df_cut_cut_cut = df_cut_cut.Filter(filter_string);
      auto h2d_var3_var4_sub = df_cut_cut_cut.Histo2D({TString(Form("h2d_%s_%s_%zu_sub", var3.Data(), var4.Data(), idx)).Data(), TString(Form("%s-%s distribution", var3.Data(), var4.Data())).Data(), (int)(var3_bins.size() - 1), var3_bins.data(), (int)(var4_bins.size() - 1), var4_bins.data()}, string(var3), string(var4));
      var3_var4_subhists_arr.push_back(h2d_var3_var4_sub);
    }

    var3_var4_subhists.push_back(var3_var4_subhists_arr);
  }

  // force z-axis to be identical to main hist
  double zmin = main_hist->GetMinimum();
  double zmax = main_hist->GetMaximum();

//   for (auto& h : var1_var2_hists) {
//     h.SetMinimum(zmin);
//     h.SetMaximum(zmax);
//   }

  std::vector<TH2D> _var1_var2_hists;
  std::vector<TH2D> _var3_var4_hists;
  std::vector<std::vector<TH2D>> _var3_var4_subhists;
  
  for (auto& hist : var1_var2_hists) {
    _var1_var2_hists.push_back(*hist);
  }
  
  for (auto& hist : var3_var4_hists) {
    _var3_var4_hists.push_back(*hist);
  }

  for (const auto subhist_list : var3_var4_subhists) {
    std::vector<TH2D> _subhist_list;
    for (auto hist : subhist_list) {
        _subhist_list.push_back(*hist);
    }
    _var3_var4_subhists.push_back(_subhist_list);
  }
  return std::make_tuple(*main_hist, _var1_var2_hists, _var3_var4_hists, _var3_var4_subhists);
}


void prepare_hist(TH2D& hist, const std::string& title, int num_contours) {
    hist.SetTitle(title.c_str());
    hist.SetContour(num_contours);
    hist.SetStats(0);
}

void gen_canvas(int boundary_idx, TH2D h_xQ2, std::vector<TH2D>& h_xQ2_list, std::vector<TH2D>& h_zpT_list, std::vector<std::vector<TH2D>>& h_zpT_sublist) {
    // Create a canvas and divide it into two pads for x-Q2 and z-pT histograms
    TCanvas c("c", "x-Q2 distribution", 1200, 600);
    c.Divide(2, 1);

    // Prepare x-Q2 histogram for selected region with log scale axes
    c.cd(1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetRightMargin(0.15);
    TH2D& h2d_xQ2 = h_xQ2_list[boundary_idx];
    prepare_hist(h2d_xQ2, ";x;Q^{2} [GeV]", 50);
    h2d_xQ2.Draw("colz");

    // Create z-pT histogram for selected region
    TH2D& h2d_zpT = h_zpT_list[boundary_idx];
    std::vector<TH2D>& h2d_zpT_sublist = h_zpT_sublist[boundary_idx];
    
    // Prepare z-pT histogram with log scale axes
    c.cd(2);
    gPad->SetRightMargin(0.15);
    prepare_hist(h2d_zpT, "z pT distribution;z;p_{T} [GeV]", 50);
    h2d_zpT.Draw("colz");

    // Draw the x-Q2 histogram with all boundaries shown
    TCanvas c2("c2", "c2", 600, 600);
    gPad->SetLogx();
    gPad->SetLogy();

    TH2D h2d_xQ2_all = h_xQ2;
    prepare_hist(h2d_xQ2_all, ";x;Q^{2} [GeV]", 50);
    h2d_xQ2_all.Draw("colz");

    c.SaveAs("c1.pdf");
    c2.SaveAs("c2.pdf");
}


int mcplotter(){
    std::vector<std::vector<std::string>> xQ2_boundaries = generate_input_boundaries({{0.07,0.12},{0.12,0.2},{0.2,0.275},{0.275,0.42},{0.42,1},
                                                                                      {0.12,0.15},{0.15,0.22},{0.22,0.29},{0.29,0.42}},
                                                                                     {0,-1,-1,-1,0,1,1,1,1});

    std::vector<std::vector<std::string>> zpT_boundaries = generate_z_pt_boundaries({0.15,0.3,0.38,0.5,0.6,0.7,0.95},{0,0.2,0.4,0.6,0.8,1,1.5});

    // Example usage
    std::vector<double> x_bins = logspace_bins(5e-2, 1, 100);
    std::vector<double> Q2_bins = logspace_bins(1, 20, 100);
    std::vector<double> z_bins = linspace(0, 1, 100);
    std::vector<double> pt_bins = linspace(-0.3, 2, 100);

    // Adjust the input file paths as needed
    //std::string infile = "../../projects/ana_v2/volatile/data/piplus_piminus/MC_RGA_3051_0.root";
    std::string infile = "../../macros/hipoLUND2tree.root";
    auto results = get_hists(infile, xQ2_boundaries, zpT_boundaries, "x", "Q2", x_bins, Q2_bins, "z", z_bins, "pTtot", pt_bins);
    TH2D xQ2_hist_main = std::get<0>(results);
    std::vector<TH2D> xQ2_hists = std::get<1>(results);
    std::vector<TH2D> zpT_hists = std::get<2>(results);
    std::vector<std::vector<TH2D>> zpT_subhists = std::get<3>(results);

    gen_canvas(0, xQ2_hist_main, xQ2_hists, zpT_hists, zpT_subhists);

    return 0;
}