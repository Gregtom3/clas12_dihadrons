from binning import *
from smearing_matrix import *
import copy
from tqdm import tqdm
ROOT.gSystem.Load("/work/clas12/users/gmat/packages/RooUnfold_cvmfs/RooUnfold/libRooUnfold.so")
#ROOT.gSystem.Load("/home/gmat/link_to_clas12/packages/RooUnfold/libRooUnfold.so")

class Unfold:
    
    def perform_RooUnfold(usable_bins,true_bins,reco_bins,data_bins,method="bayes"):
        # Load the C++ function into ROOT
        ROOT.gROOT.ProcessLine('.L ./unfold.C')


        # Create a ROOT std::vector<int> object
        tb_vector = ROOT.std.vector('int')()
        rb_vector = ROOT.std.vector('int')()
        db_vector = ROOT.std.vector('int')()
        # Fill the vector with the elements of the Python list
        print("LOADING MONTE CARLO VECTORS")
        for i in tqdm(range(len(true_bins))):
            tb_vector.push_back(int(true_bins[i]))
            rb_vector.push_back(int(reco_bins[i]))
        print("LOADING EXPERIMENTAL DATA VECTORS")
        for j in tqdm(range(len(data_bins))):
            db_vector.push_back(int(data_bins[j]))
        # Perform the unfolding
        print("PERFORMING UNFOLDING")
        ROOT.perform_RooUnfold(usable_bins, tb_vector, rb_vector, db_vector, "bayes")
        # Open the saved unfolded hists file and return the histograms
        f = ROOT.TFile.Open("/work/clas12/users/gmat/tmp/unfolded_hists.root")
        reco_hist = f.Get("reco_hist")
        true_hist = f.Get("true_hist")
        data_hist = f.Get("data_hist")
        unfolded_hist = f.Get("unfolded_hist")
        return copy.deepcopy(reco_hist),copy.deepcopy(true_hist),copy.deepcopy(data_hist),copy.deepcopy(unfolded_hist)
        
#     def perform_RooUnfold(usable_bins,true_bins,reco_bins,data_bins,method="bayes"):
#         # Unfold everything except under/overflow
        
#         bw = 0.5 # bin-width
        
#         response = ROOT.RooUnfoldResponse (usable_bins,1.0,usable_bins+1.0)
#         true_hist = ROOT.TH1D("true_hist","",usable_bins,1.0,usable_bins+1.0)
#         reco_hist = ROOT.TH1D("reco_hist","",usable_bins,1.0,usable_bins+1.0)
#         data_hist = ROOT.TH1D("data_hist","",usable_bins,1.0,usable_bins+1.0)

        
#         for tb,rb in zip(true_bins,reco_bins):
#             if(tb==0 or rb==0): # Skip passed the underflow/overflow bins
#                 continue
#             if(tb==-1 and rb==-1): # Skip event
#                 continue
#             if(tb!=-1): # Generated a dihadron
#                 true_hist.Fill(tb+bw)
#             if(rb==-1 and tb!=-1): # Generated a dihadron, yet it was not reconstructed
#                 response.Miss(tb+bw)
#             elif(rb!=-1 and tb==-1): # Reconstructed a dihadron where none was generated
#                 response.Fake(rb+bw)
#                 reco_hist.Fill(rb+bw)
#             else:
#                 response.Fill(rb+bw,tb+bw)
#                 reco_hist.Fill(rb+bw)
        
#         for db in data_bins:
#             if(db==0):
#                 continue
#             data_hist.Fill(db+bw)
            
#         if method == "bayes":
#             uf = ROOT.RooUnfoldBayes(response, data_hist, 25);    #  OR\
#         elif method == "svd":
#             uf= ROOT.RooUnfoldSvd (response, data_hist,20);
            
#         uf.HandleFakes(1)
#         hUnfold = uf.Hunfold();
        
#         return copy.deepcopy(reco_hist),copy.deepcopy(true_hist),copy.deepcopy(data_hist),copy.deepcopy(hUnfold)
    
    def plot_unfolded(reco_hist=None, true_hist=None, unfolded_hist=None):
        canvas = ROOT.TCanvas("canvas", "Histograms", 800, 600)
        ROOT.gPad.SetLeftMargin(0.15)

        legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)

        if reco_hist is not None:
            reco_clone = reco_hist.Clone()
            reco_clone.SetLineColor(ROOT.kRed)
            reco_clone.SetTitle(";Bin Entry;Counts")
            reco_clone.Draw("hist same")
            legend.AddEntry(reco_clone, "Reco Hist", "l")

        if true_hist is not None:
            true_clone = true_hist.Clone()
            true_clone.SetLineColor(ROOT.kBlue)
            true_clone.Draw("hist same")
            legend.AddEntry(true_clone, "True Hist", "l")

        if unfolded_hist is not None:
            unfold_clone = unfolded_hist.Clone()
            unfold_clone.SetLineColor(ROOT.kGreen)
            unfold_clone.Draw("hist same")
            legend.AddEntry(unfold_clone, "Unfold Hist", "l")

        legend.Draw()
        canvas.Draw()

        return copy.deepcopy(canvas)
    
    
    
    def plot_pulls(reco_hist, true_hist, unfolded_hist):
        # Create a canvas to display the plots
        canvas = ROOT.TCanvas("canvas", "Pulls", 800, 600)
        ROOT.gStyle.SetOptStat(0)
        # Calculate the pulls by subtracting true_hist from reco_hist and unfolded_hist
        pulls_reco = reco_hist.Clone()
        pulls_reco.Add(true_hist, -1)
        pulls_reco.SetName("pulls_reco")
        pulls_reco.SetTitle("Pulls (reco)")
        pulls_reco.SetLineColor(ROOT.kRed)
        pulls_reco.SetMarkerColor(ROOT.kRed)
        pulls_reco.SetMarkerStyle(ROOT.kFullCircle)

        pulls_unfold = unfolded_hist.Clone()
        pulls_unfold.Add(true_hist, -1)
        pulls_unfold.SetName("pulls_unfold")
        pulls_unfold.SetTitle("Pulls (unfold)")
        pulls_unfold.SetLineColor(ROOT.kBlue)
        pulls_unfold.SetMarkerColor(ROOT.kBlue)
        pulls_unfold.SetMarkerStyle(ROOT.kFullSquare)

        # Set the bin content errors as the square root of bin counts for all histograms
        for bin in range(1, true_hist.GetNbinsX() + 1):
            reco_content = reco_hist.GetBinContent(bin)
            reco_error = reco_content ** 0.5
            reco_hist.SetBinError(bin, reco_error)

            true_content = true_hist.GetBinContent(bin)
            true_error = true_content ** 0.5
            true_hist.SetBinError(bin, true_error)

            unfold_content = unfolded_hist.GetBinContent(bin)
            unfold_error = unfold_content ** 0.5
            unfolded_hist.SetBinError(bin, unfold_error)

            # Calculate the pulls with error propagation (avoiding division by zero)
            if true_content != 0:
                pull_reco = pulls_reco.GetBinContent(bin)
                pull_reco_error = ((reco_error / true_content) ** 2 + (reco_content * true_error / true_content ** 2) ** 2) ** 0.5
                pulls_reco.SetBinError(bin, pull_reco_error)

                pull_unfold = pulls_unfold.GetBinContent(bin)
                pull_unfold_error = ((unfold_error / true_content) ** 2 + (unfold_content * true_error / true_content ** 2) ** 2) ** 0.5
                pulls_unfold.SetBinError(bin, pull_unfold_error)
            else:
                pulls_reco.SetBinContent(bin, 0)
                pulls_reco.SetBinError(bin, 0)
                pulls_unfold.SetBinContent(bin, 0)
                pulls_unfold.SetBinError(bin, 0)

        # Divide the pulls histograms by true_hist for error propagation
        pulls_reco.Divide(true_hist)
        pulls_unfold.Divide(true_hist)

        # Draw the pull histograms on the canvas, excluding data points with error bars larger than epsilon
        epsilon = 0.4
        for bin in range(1, true_hist.GetNbinsX() + 1):
            pull_reco_error = pulls_reco.GetBinError(bin)
            if pull_reco_error >= epsilon:
                pulls_reco.SetBinContent(bin, 0)
                pulls_reco.SetBinError(bin, 0)

            pull_unfold_error = pulls_unfold.GetBinError(bin)
            if pull_unfold_error >= epsilon:
                pulls_unfold.SetBinContent(bin, 0)
                pulls_unfold.SetBinError(bin, 0)

        pulls_reco.Draw("PE")
        pulls_unfold.Draw("PE SAME")

        pulls_reco.GetYaxis().SetRangeUser(-1, 1)
        # Create a legend
        legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
        legend.AddEntry(pulls_reco, "Pulls (reco)", "lep")
        legend.AddEntry(pulls_unfold, "Pulls (unfold)", "lep")
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.Draw()

        # Update the canvas and display the plot
        canvas.Update()
        canvas.Draw()

        return copy.deepcopy(canvas)