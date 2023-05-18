from binning import *
from smearing_matrix import *
import copy

ROOT.gSystem.Load("/home/gmat/link_to_clas12/packages/RooUnfold/libRooUnfold.so")

class Unfold:
        
    def perform_RooUnfold(smearing_matrix,method):
        
        usable_bins = smearing_matrix.bin_manager.total_bins - 1 # Unfold everything except under/overflow
        
        response = ROOT.RooUnfoldResponse (usable_bins,1.0,usable_bins+1.0)
        
        true_bins = smearing_matrix.true_bins
        reco_bins = smearing_matrix.reco_bins
        
        for tb,rb in zip(true_bins,reco_bins):
            response.Fill(rb,tb)
            
        true_vector = np.sum(smearing_matrix.smearing_matrix, axis=1)
        reco_vector = np.sum(smearing_matrix.smearing_matrix, axis=0)
        
        true_hist = ROOT.TH1D("true_hist","",usable_bins,1.0,usable_bins+1.0)
        reco_hist = ROOT.TH1D("reco_hist","",usable_bins,1.0,usable_bins+1.0)
        
        for i,(x,xt) in enumerate(zip(reco_vector,true_vector)):
            if i==0:
                continue
            else:
                reco_hist.Fill(i,x)
                true_hist.Fill(i,xt)
        
        if method == "bayes":
            uf = ROOT.RooUnfoldBayes(response, reco_hist, 4);    #  OR\
        elif method == "svd":
            uf= ROOT.RooUnfoldSvd (response, reco_hist,20);
            
        hUnfold = uf.Hunfold();
        
        return copy.deepcopy(reco_hist),copy.deepcopy(true_hist),copy.deepcopy(hUnfold)