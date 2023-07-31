import uproot
import numpy as np
from array import array
import re
import ROOT
import matplotlib.pyplot as plt
import yaml
import os
import copy
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Pool
ROOT.gErrorIgnoreLevel = ROOT.kFatal
from natsort import natsorted # pip install natsort
from tqdm import tqdm





class BinManager:
    def __init__(self,u):
        self.u = u
        self.bins = []
        self.bin_names = []
        #print("Loading MCmatch")
        self.MCmatch = self.u["MCmatch"].array(library="np")

    def create_bin(self,bin_func):
        pattern = ":(.*?):"
        self.bin_names += re.findall(pattern, bin_func)
        self.bin_names = list(set(self.bin_names)) # Remove duplicates
        eval_bin_func = bin_func # To be edited in the for loop
        for bin_name in self.bin_names:
            try:
                exec(f"self.b_{bin_name}")
            except:
                #print(f"Creating bin variable `self.b_{bin_name}`")
                exec(f"self.b_{bin_name}=self.u[\"{bin_name}\"].array(library=\"np\")")
            eval_bin_func = eval_bin_func.replace(f":{bin_name}:",f"self.b_{bin_name}")
            
            
        include_indecies = eval(eval_bin_func) * (self.MCmatch==1)
        self.bins.append(Bin(bin_func,include_indecies))

    def get_bin_avg_of_branch(self,branch_values):
        return np.mean(branch_values[self.include_indecies])
    
    def get_bin_avgs(self):
        avgs=[]
        for bin_name in self.bin_names:
            avgs.append({"name": bin_name,
                         "x"   : []})
            d = avgs[-1]
            for myBin in self.bins:
                exec(f"d['x'].append(np.mean(self.b_{bin_name}[myBin.include_indecies]))")
        return avgs
    
class Bin(BinManager):
    def __init__(self,bin_func,include_indecies):
        self.bin_func = bin_func
        self.neat_bin_func = self.parse_bin_func(bin_func)
        self.include_indecies = include_indecies
        self.sideband =  {"var": "",
                          "purity": 0,
                          "signal":     {"vmin": [],
                                         "vmax": [],
                                         "fitpars": [],
                                         "fiterrs": []},
                          "signal+bkg": {"vmin": [],
                                         "vmax": [],
                                         "fitpars": [],
                                         "fiterrs": []},
                          "bkg":        {"vmin": [],
                                         "vmax": [],
                                         "fitpars": [],
                                         "fiterrs": []}}
        # The below two lists are populated when sidebands are not used
        self.fitpars  = []
        self.fiterrs  = [] 
    
    @staticmethod
    def parse_bin_func(bf):
        pattern = r':(\w+):'
        return re.sub(pattern, r'\1', bf)
    
    def print(self):
        
        # Print the Bin Function
        print(f"Bin ({self.neat_bin_func}):",("" if self.sideband["purity"]==0 else f"Purity = {self.sideband['purity']}"))
        
        # Code for printing fit variables
        def print_vars(vs,es,sp):
            if sp != "":
                print(f"\t Sideband Region:")
            for i,[v,e] in enumerate(zip(vs,es)):
                
                print(f"\t\t[{i}]({sp}) = {v:.4f} +/- {e:.4f}")
        
        # Code for printing sideband/region specific asymmetries
        if self.sideband["var"]!="": # Print out for each region
            for species in ["signal","signal+bkg","bkg"]:
                
                vmin = self.sideband[species]["vmin"]
                vmax = self.sideband[species]["vmax"]

                fitpars = self.sideband[species]["fitpars"]
                fiterrs = self.sideband[species]["fiterrs"]

                print_vars(fitpars,fiterrs,species)

        else:
            print_vars(self.fitpars,self.fiterrs,"")
            
            
class Injector:
    ####################################################################################
    # Intializer + Initializer Functions
    ####################################################################################
    def __init__(self, pion_pair, filename, outdir):
        self.filename = filename
        self.u = uproot.open(self.filename)
        self.u = self.u["dihadron_cuts"]
        self.N = self.u.num_entries
        self.pion_pair = pion_pair
        self.outdir = outdir
        self.set_pion_pids()
        self.load_azi_vars()
        self.weight_funcs = [ ["0.2","0.0","sin(:phi_R0:)"] ]
        self.bPol = 1.0
        self.nbins = 30 # For 2d fits
        self.sideband =  {"var": "",
                          "purity": 0,
                          "signal":     {"vmin": [],
                                         "vmax": [],
                                         "fitpars": [],
                                         "fiterrs": []},
                          "signal+bkg": {"vmin": [],
                                         "vmax": [],
                                         "fitpars": [],
                                         "fiterrs": []},
                          "bkg":        {"vmin": [],
                                         "vmax": [],
                                         "fitpars": [],
                                         "fiterrs": []}}
        self.bin_manager = BinManager(self.u)
        self.set_true_dihadron_indecies()

    # Load commonly used azimuthal variables
    def load_azi_vars(self):
        #print("Loading PhiH")
        self.PhiH = self.u["phi_h"].array(library="np")
        self.TruePhiH = self.u["truephi_h"].array(library="np")
        #print("Loading PhiR0")
        self.PhiR0 = self.u["phi_R0"].array(library="np")
        self.TruePhiR0 = self.u["truephi_R0"].array(library="np")
        #print("Loading PhiR1")
        self.PhiR1 = self.u["phi_R1"].array(library="np")
        self.TruePhiR1 = self.u["truephi_R1"].array(library="np")

    # Determine the corresponding pids for the pion pairs
    def set_pion_pids(self):
        if self.pion_pair == "piplus_piplus":
            self.pids = [211,211]
        elif self.pion_pair == "piplus_pi0":
            self.pids = [211,111]
        elif self.pion_pair == "piplus_piminus":
            self.pids = [211,-211]
        elif self.pion_pair == "pi0_pi0":
            self.pids = [111, 111]
        elif self.pion_pair == "piminus_pi0":
            self.pids = [-211, 111]
        elif self.pion_pair == "piminus_piminus":
            self.pids = [-211, -211]
        else:
            # Handle an unknown or unsupported pion pair
            raise ValueError("Unsupported pion pair: {}".format(self.pion_pair))
    
    # Determine the Monte Carlo indecies containing a true dihadron match
    def set_true_dihadron_indecies(self):
        #print("Loading indecies containing true dihadrons")
        truepid_1 = self.u["truepid_1"].array(library="np")
        truepid_2 = self.u["truepid_2"].array(library="np")
        trueparentpid_1 = self.u["trueparentpid_1"].array(library="np")
        trueparentpid_2 = self.u["trueparentpid_2"].array(library="np")
        
        true_dihadron_indecies = np.ones(self.N)
        
        if(self.pids[0]==111):
            true_dihadron_indecies *= (trueparentpid_1==111)
        else:
            true_dihadron_indecies *= (truepid_1==self.pids[0])
            
        if(self.pids[1]==111):
            true_dihadron_indecies *= (trueparentpid_2==111)
        else:
            true_dihadron_indecies *= (truepid_2==self.pids[1])
            
        self.true_dihadron_indecies = (true_dihadron_indecies==1)
        self.false_dihadron_indecies = (true_dihadron_indecies==0)
        
    ####################################################################################
    # String parsing Functions
    ####################################################################################
    @staticmethod
    def parse_func(s):
        pattern = r':(\w+):'
        new_s = s
        new_s = new_s.replace(":phi_h:","self.PhiH").replace(":phi_R0:","self.PhiR0").replace(":phi_R1:","self.PhiR1")
        new_s = re.sub(pattern, r'self.u["\1"].array(library="np")', new_s)
        new_s = new_s.replace("sin","np.sin")
        return new_s
    
    @staticmethod
    def parse_true_func(s):
        pattern = r':(\w+):'
        new_s = s
        new_s = new_s.replace(":phi_h:","self.TruePhiH").replace(":phi_R0:","self.TruePhiR0").replace(":phi_R1:","self.TruePhiR1")
        new_s = re.sub(pattern, r'self.u["true\1"].array(library="np")', new_s)
        new_s = new_s.replace("sin","np.sin")
        return new_s
    
    @staticmethod
    def parse_func_latex(s):
        new_s = s.replace(":phi_h:","\phi_{h}").replace(":phi_R0:","\phi_{RT}").replace(":phi_R1:","\phi_{R\perp}")
        return new_s
    
    ####################################################################################
    # User-Input Functions
    ####################################################################################
    # Load in asymmetry injection
    # We feed it a list of lists
    # The inner lists contain 3 elements, the signal injection, the background injection, and the modulation
    def load_weight_funcs(self, weight_funcs):
        self.weight_funcs = weight_funcs
    
    # Set the beam polarization [0,1]
    def set_beam_polarization(self, bPol):
        self.bPol = bPol
        
    # Set the variable used to allocate sideband regions
    # This is typically "M2" for the pi0 analyses
    def set_sideband_var(self, var):
        self.sideband["var"]=var
        #print("Loading sideband variable",var)
        self.sideband_var = self.u[var].array(library="np")
        
    # Create a sideband region
    # species --> name of region
    # vmin,vmax --> region boundaries
    def addSideband(self,species,vmin,vmax):
        if(species not in self.sideband.keys() or species=="signal"):
            raise ValueError("Species",species,"not known...\n Options are [signal+bkg] [bkg]")
        
        self.sideband[species]["vmin"].append(vmin)
        self.sideband[species]["vmax"].append(vmax)
        
        if species=="signal+bkg":
            self.sideband["signal"]["vmin"].append(vmin)
            self.sideband["signal"]["vmax"].append(vmax)
        
    # Create a new bin
    # ex: bin_func = "(:Mh: > 0.03) & (:Mh: < 0.06)"
    # Colons denote variable names
    # NOTE: Parentheses are important here
    def addBin(self,bin_func):
        self.bin_manager.create_bin(bin_func)
    
    def create_bins_linspace(self, binname, start, end, step):
        bin_start = start
        while bin_start < end:
            bin_end = min(bin_start + step, end)
            bin_range = f"(:{binname}: > {bin_start}) & (:{binname}: < {bin_end})"
            self.addBin(bin_range)
            bin_start = bin_end
            
    def create_bins(self, binname, bin_edges):
        for i in range(len(bin_edges) - 1):
            bin_start = bin_edges[i]
            bin_end = bin_edges[i+1]
            bin_range = f"(:{binname}: > {bin_start}) & (:{binname}: < {bin_end})"
            self.addBin(bin_range)
            
    ####################################################################################
    # Internal functions
    ####################################################################################
    # Define the helicities per the user-inputted asymmetries
    def get_helicities(self):
        
        weight_funcs = self.weight_funcs
        
        def parse_weight_func():
            full_sig_func = "0.5+0.5*("
            full_bkg_func = "0.5+0.5*("
            for weight_func in self.weight_funcs:
                asym_sig = weight_func[0]; asym_sig=self.parse_true_func(asym_sig)
                asym_bkg = weight_func[1]; asym_bkg=self.parse_true_func(asym_bkg)
                mod  = weight_func[2]; mod=self.parse_true_func(mod)
                full_sig_func += "+(" + asym_sig + ")*(" + mod + ")"
                full_bkg_func += "+(" + asym_bkg + ")*(" + mod + ")"
            full_sig_func += ")"
            full_bkg_func += ")"
            return full_sig_func, full_bkg_func
        
        def get_hel():
            # Generate a list of random numbers between 0 and 1
            hel = np.random.rand(self.N)
            
            # Store the 50/50 helicity indices
            idx_5050 = np.random.rand(self.N)>self.bPol
            
            # Determine the probability weighting string
            new_weight_sig_func, new_weight_bkg_func = parse_weight_func()
            
            # Assign helicity to -1 or 1 depending on probability
            # Also consider if this dihadron is signal or background
            hel[self.true_dihadron_indecies] = eval("2*(hel<"+new_weight_sig_func+")-1")[self.true_dihadron_indecies]
            hel[self.false_dihadron_indecies] = eval("2*(hel<"+new_weight_bkg_func+")-1")[self.false_dihadron_indecies]
            
            # Create a random list of -1 or +1 helicities
            rand_hel = np.random.choice([-1,1],self.N) 
            
            # For the 50/50 indices, assign them to -1 or +1 at random
            hel[idx_5050] = rand_hel[idx_5050]
            
            return hel
    
        self.hel  = get_hel()
        
        
    
    # Draw helicity dependent histograms to illustrate injection
    def draw_hel_hists(self):
        
        Nplots = len(self.weight_funcs)+3
        Ncols = 3
        Nrows = int(np.ceil(Nplots/Ncols))
        
        fig,axs = plt.subplots(Nrows,Ncols,dpi=150,figsize=(3*Ncols,3*Nrows))
        axs=axs.flatten()
        r,c = 0,0
        for i,wf in enumerate(self.weight_funcs):
            
            func = wf[2]
            
            x = eval(self.parse_func(func))

            # Draw +1 helicity
            axs[i].hist(x[self.hel==1],bins=100,histtype="step",color="red",label="Helicity = +1")

            # Draw -1 helicity
            axs[i].hist(x[self.hel==-1],bins=100,histtype="step",color="blue",label="Helicity = -1")

            axs[i].set_xlabel("$"+self.parse_func_latex(func)+"$")

            axs[i].legend()
        
        def style_2d(ax,h=""):
            ax.set_xlabel("$\phi_{h}$")
            ax.set_ylabel("$\phi_{R}$")
            text_props = dict(
                horizontalalignment='center',
                verticalalignment='center',
                bbox=dict(
                    facecolor='white',
                    edgecolor='gray',
                    boxstyle='round'
                )
            )

            # Add the text with white background and gray border
            if h:
                ax.text(0.75, 0.88, f"Helicity={h}", transform=ax.transAxes, **text_props)
        
        ##################################
        # Create the phiR vs phiH hists
        ##################################
        axs[i+1].hist2d(self.PhiH[self.hel==1],self.PhiR0[self.hel==1],bins=self.nbins,range=((-np.pi,np.pi),(-np.pi,np.pi)), cmap='rainbow')
        style_2d(axs[i+1],"+1")
        axs[i+2].hist2d(self.PhiH[self.hel==-1],self.PhiR0[self.hel==-1],bins=self.nbins,range=((-np.pi,np.pi),(-np.pi,np.pi)), cmap='rainbow')
        style_2d(axs[i+2],"-1")
        
        ##################################
        # Create the asym phiR vs phiH hist
        ##################################
        hist_1, xedges, yedges = np.histogram2d(self.PhiH[self.hel==1],self.PhiR0[self.hel==1],bins=self.nbins,range=((-np.pi,np.pi),(-np.pi,np.pi)))
        hist_2, _, _ = np.histogram2d(self.PhiH[self.hel==-1],self.PhiR0[self.hel==-1],bins=self.nbins,range=((-np.pi,np.pi),(-np.pi,np.pi)))

        # Calc asym
        hist_diff = (hist_1 - hist_2)/(hist_1 + hist_2)

        # Generate a grid for the plot
        X, Y = np.meshgrid(xedges, yedges)

        # Plot the difference using pcolormesh or imshow
        c = axs[i+3].pcolormesh(X, Y, hist_diff.T, cmap='rainbow') 
        style_2d(axs[i+3])
        plt.colorbar(c, ax=axs[i+3])
        
        
        plt.tight_layout()

        # Remove empty subplots
        for j in range(Nrows*Ncols-1, i+3,-1):
                fig.delaxes(axs[j])
        
    
    # Calculate the asymmetries depending on the "include_indecies"
    # By setting include indecies to be all 1, all events are used
    # when calculating the asymmetry
    def calculate_asymmetry(self,include_indecies):
    
        # Initialize the helicity +1 and helicity -1 histograms
        nbins = 30

        HPOS = ROOT.TH2F("HPOS","",self.nbins,-np.pi,np.pi,self.nbins,-np.pi,np.pi)
        HNEG = ROOT.TH2F("HNEG","",self.nbins,-np.pi,np.pi,self.nbins,-np.pi,np.pi)

        # Create the TF2
        fit_string = ""

        for i,wf in enumerate(self.weight_funcs):
            s = wf[2]
            s = s.replace(":phi_h:","x").replace(":phi_R0:","y").replace(":phi_R1:","y")
            fit_string+=f"+[{i}]*{s}"

        tf2 = ROOT.TF2("tf2",fit_string,-np.pi,np.pi,-np.pi,np.pi)

        # Fill the histograms
        #for phiH,phiR,h,inc in zip(self.PhiH,self.PhiR0,self.hel,include_indecies):
        for phiH,phiR,h,inc in zip(self.TruePhiH,self.TruePhiR0,self.hel,include_indecies):
            if inc != 1: # Skip undesired region
                continue
            if h == 1:
                HPOS.Fill(phiH,phiR)
            else:
                HNEG.Fill(phiH,phiR)
        
        # Calculate asym histogram
        HNUM = HPOS.Clone()
        HDEN = HPOS.Clone()
        
        HNUM.Sumw2(True)
        HDEN.Sumw2(True)
        
        HNUM.Add(HNEG,-1)
        HDEN.Add(HNEG,+1)
        
        H=HNUM.Clone()
        H.Divide(HDEN)
        
        H.Scale(1/self.bPol)
        
        # Fit the histogram
        H.Fit(tf2,"NQR0")
    
        # Return the fit parameters and fit errors
        return np.array([tf2.GetParameter(i) for i in range(len(self.weight_funcs))]), np.array([tf2.GetParError(i) for i in range(len(self.weight_funcs))])
    
    # For each user defined 
    def calculate_asymmetries_for_each_sideband(self,myBin):
        
        myBin.sideband["purity"] = self.calculate_purity(myBin)
        
        for species in ["signal+bkg","bkg"]:

            # Find the indecies where we should calculate the asymmetries
            vmin = myBin.sideband[species]["vmin"]
            vmax = myBin.sideband[species]["vmax"]
            
            # Initialize an empty boolean array of the same length as self.sideband_var, all False
            include_indecies = np.full(len(self.sideband_var), False)

            # Assume vmin and vmax are lists of same length
            for vmin_i, vmax_i in zip(vmin, vmax):
                include_indecies |= (myBin.include_indecies) & (self.sideband_var>vmin_i) & (self.sideband_var<vmax_i)
            
            pars,errs = self.calculate_asymmetry(include_indecies)

            myBin.sideband[species]["fitpars"] = pars
            myBin.sideband[species]["fiterrs"] = errs
        
        # Now we calculate the full signal pars
        u = myBin.sideband["purity"]
        SB = myBin.sideband["signal+bkg"]["fitpars"]
        eSB = myBin.sideband["signal+bkg"]["fiterrs"]
        B = myBin.sideband["bkg"]["fitpars"]
        eB = myBin.sideband["bkg"]["fiterrs"]
        
        myBin.sideband["signal"]["fitpars"] = (1/u) * SB - ((1-u)/(u))*B
        myBin.sideband["signal"]["fiterrs"] = np.sqrt((1/u)**2 * eSB**2 + ((1-u)/(u))**2*eB**2)
        

    # Print out fit details
    def print(self):
        print("="*15,"Fit Results","="*15)
        for myBin in self.bin_manager.bins:
            myBin.print()
            print("\n")
    
    # Calculate the purity of region
    def calculate_purity(self,myBin):
        
        c=ROOT.TCanvas()
        
        include_indecies = myBin.include_indecies
        
        H = ROOT.TH1F("H",myBin.sideband["var"],100,0.03,0.4)
        
        for ind,v in zip(include_indecies,self.sideband_var):
            if ind==0:
                continue
            else:
                H.Fill(v)
                
        
        max_bin_center = -np.inf
        for i in range(1,H.GetNbinsX()+1):
            if H.GetBinContent(i)>10:
                max_bin_center = H.GetBinCenter(i)
        
        H.Scale(1/H.Integral())

            
        tf1 = ROOT.TF1("tf1","gaus(0)+pol4(3)",0.03,max_bin_center)
        
        tf1.SetParLimits(1,0.125,0.14)
        tf1.SetParLimits(2,0.001,0.02)
        
        H.Fit(tf1,"NQR0")
        
        tf1_1 = ROOT.TF1("tf1_1","gaus(0)",0.03,0.4)
        tf1_2 = ROOT.TF1("tf1_2","pol4(0)",0.03,0.4)
        for i in range(3):
            tf1_1.SetParameter(i,tf1.GetParameter(i))
        for i in range(3,9):
            tf1_2.SetParameter(i-3,tf1.GetParameter(i))
        Nsig=tf1_1.Integral(0.106,0.166)
        Nbkg=tf1_2.Integral(0.106,0.166)
        
        
        H.Draw("hist")
        tf1.Draw("same")
        c.SaveAs(f"{self.outdir}/canv___{myBin.neat_bin_func}.png")
        
        return Nsig/(Nsig+Nbkg)
    
    # Plot the results as a function of the binning
    def plot_asym(self):
        # Get a dictionary of the bin averages
        data = self.bin_manager.get_bin_avgs()[0]
        
        Ncol = len(self.bin_manager.bins[0].fitpars)

        fig,axs = plt.subplots(1,Ncol,dpi=200,figsize=(3*Ncol,4),sharex=True,sharey=True)
        
        axs=axs.flatten()

        for col in range(Ncol):
            weight_func = self.weight_funcs[col]
            asym_sig = weight_func[0]; asym_sig=self.parse_true_func(asym_sig)
            asym_sig_values = np.array(eval(asym_sig))
            asym_bkg = weight_func[1]; asym_bkg=self.parse_true_func(asym_bkg)
            asym_bkg_values = np.array(eval(asym_bkg))
            
            if type(asym_sig_values)!=list:
                asym_sig_values = asym_sig_values * np.ones(self.N)
            if type(asym_bkg_values)!=list:
                asym_bkg_values = asym_bkg_values * np.ones(self.N)
            
            ax = axs[col]

            ax.set_xlabel(data["name"])

            x = data["x"]
            y = [myBin.fitpars[col] for myBin in self.bin_manager.bins]
            yerr = [myBin.fiterrs[col] for myBin in self.bin_manager.bins]

            ax.errorbar(x,y,yerr=yerr,capsize=3,fmt="ko",label="Measured")
            ax.grid()
            ax.axhline(0,color="black",linestyle="dashed")

            ax.set_ylabel(f"$A_{{LU}}^{{{self.parse_func_latex(weight_func[2])}}}$",fontsize=15)
            trueysig = [myBin.get_bin_avg_of_branch(asym_sig_values) for myBin in self.bin_manager.bins]
            ax.plot(x,trueysig,"r*",markersize=10,label="Injected (sig)")
            trueybkg = [myBin.get_bin_avg_of_branch(asym_bkg_values) for myBin in self.bin_manager.bins]
            ax.plot(x, trueybkg, marker="*", markerfacecolor="white", markeredgecolor="red", linestyle="", markersize=10, label="Injected (bkg)")
            ax.legend()
            
                
        ymax = np.amax(np.abs(ax.get_ylim()))
        if ymax>1:
            ymax = 1
        ax.set_ylim(-ymax,ymax)
        
        plt.tight_layout()
        
    # Get the asym_sig_values and asym_bkg_values for a given row
    def get_asym_values(self,row):
        weight_func = self.weight_funcs[row]
        asym_sig = weight_func[0]; asym_sig=self.parse_true_func(asym_sig)
        asym_sig_values = np.array(eval(asym_sig))
        asym_bkg = weight_func[1]; asym_bkg=self.parse_true_func(asym_bkg)
        asym_bkg_values = np.array(eval(asym_bkg))

        if type(asym_sig_values)!=list:
            asym_sig_values = asym_sig_values * np.ones(self.N)
        if type(asym_bkg_values)!=list:
            asym_bkg_values = asym_bkg_values * np.ones(self.N)
            
        return asym_sig_values, asym_bkg_values
    # Plot the results as a function of the binning
    def plot_asym_3(self):
        # Get a dictionary of the bin averages
        data = self.bin_manager.get_bin_avgs()[0]
        
        Nrow = len(self.bin_manager.bins[0].sideband["signal"]["fitpars"])

        fig,axs = plt.subplots(Nrow,3,dpi=200,figsize=(12,3*Nrow),sharex=True,sharey=True)
        plt.subplots_adjust(wspace=0.6)
        
        axs=axs.flatten()

        for row in range(Nrow):
            asym_sig_values, asym_bkg_values = get_asym_values(row)
            
            if type(asym_sig_values)!=list:
                asym_sig_values = asym_sig_values * np.ones(self.N)
            if type(asym_bkg_values)!=list:
                asym_bkg_values = asym_bkg_values * np.ones(self.N)
                
            for ax,species in zip(axs[3*row:row*3+3],["signal","signal+bkg","bkg"]):
                ax.set_xlabel(data["name"])
                if row == 0:
                    ax.set_title(species+" BSAs")
                
                x = data["x"]
                y = [myBin.sideband[species]["fitpars"][row] for myBin in self.bin_manager.bins]
                yerr = [myBin.sideband[species]["fiterrs"][row] for myBin in self.bin_manager.bins]
                
                ax.errorbar(x,y,yerr=yerr,capsize=3,fmt="ko",label="Measured")
                ax.grid()
                ax.axhline(0,color="black",linestyle="dashed")
                if species=="signal":
                    ax.set_ylabel(f"$A_{{LU}}^{{{self.parse_func_latex(weight_func[2])}}}$",fontsize=15)
                    truey = [myBin.get_bin_avg_of_branch(asym_sig_values) for myBin in self.bin_manager.bins]
                    ax.plot(x,truey,"r*",label="Injected")
                    ax.legend()
                elif species=="bkg":
                    truey = [myBin.get_bin_avg_of_branch(asym_bkg_values) for myBin in self.bin_manager.bins]
                    ax.plot(x,truey,"r*")
                elif species=="signal+bkg":
                    uy = [myBin.sideband["purity"] for myBin in self.bin_manager.bins]
                    ax2=ax.twinx()
                    ax2.tick_params(axis='y',colors='green')
                    ax2.plot(x,uy,"go")
                    ax2.set_ylabel("Purity", rotation=-90, labelpad=13, color='green')
                    ax2.set_ylim(0,1)
                
            ymax = np.amax(np.abs(ax.get_ylim()))
            if ymax>1:
                ymax = 1
            ax.set_ylim(-ymax,ymax)
        
                
        #plt.tight_layout()
    
    def to_yaml(self,yaml_name = "out.yaml"):
        

        data = self.bin_manager.get_bin_avgs()[0] # Get bin centers (x-values)
        x = data["x"]
        outyaml = {"x":[float(X) for X in x],
                   "xname": self.bin_manager.bin_names[0]}
        for iwf,wf in enumerate(self.weight_funcs):
            outyaml[f"mod_{iwf}"] = {"modname":self.parse_func_latex(wf[2])}
            asym_sig_values, asym_bkg_values = self.get_asym_values(iwf)
            
            for iBin,myBin in enumerate(self.bin_manager.bins):
                outyaml[f"mod_{iwf}"][f"bin_{iBin}"] = {"all":{}}
                if myBin.sideband["var"]=="": # No sideband region made
                    outyaml[f"mod_{iwf}"][f"bin_{iBin}"]["all"]["y"] = myBin.fitpars.tolist()[iwf]
                    outyaml[f"mod_{iwf}"][f"bin_{iBin}"]["all"]["yerr"] = myBin.fiterrs.tolist()[iwf]
                    outyaml[f"mod_{iwf}"][f"bin_{iBin}"]["all"]["true_signal_y"] = float(myBin.get_bin_avg_of_branch(asym_sig_values))
                    outyaml[f"mod_{iwf}"][f"bin_{iBin}"]["all"]["true_bkg_y"] = float(myBin.get_bin_avg_of_branch(asym_bkg_values))

                else:
                    for species in ["signal","signal+bkg","bkg"]:
                        outyaml[f"mod_{iwf}"][f"bin_{iBin}"][species] = {}
                        outyaml[f"mod_{iwf}"][f"bin_{iBin}"][species]["y"] = myBin.sideband[species]["fitpars"].tolist()[iwf]
                        outyaml[f"mod_{iwf}"][f"bin_{iBin}"][species]["yerr"] = myBin.sideband[species]["fiterrs"].tolist()[iwf]

                        if species == "signal":
                            outyaml[f"mod_{iwf}"][f"bin_{iBin}"][species]["truey"] = float(myBin.get_bin_avg_of_branch(asym_sig_values)) 
                        elif species == "bkg":
                            outyaml[f"mod_{iwf}"][f"bin_{iBin}"][species]["truey"] = float(myBin.get_bin_avg_of_branch(asym_bkg_values))
        
        with open(yaml_name, 'w') as file:
            yaml.dump(outyaml, file, default_flow_style=False)
            
    ####################################################################################
    # Execution command
    ####################################################################################            
    # Calculate the asymmetry histograms
    def run(self, draw_hel = False, print_out = False, draw_asym = False, print_to_yaml = False, yaml_name = ""):
        
        # Produce helicities by weighting procedure
        #print("Getting helicities...")
        self.get_helicities() 
        
        # Produce histograms of asymmetry
        if draw_hel:
            print("Drawing figures...")
            self.draw_hel_hists()
        
        # Calculate the asymmetry for each bin
        # If no bins were created, make a dummy placeholder
        if self.bin_manager.bins == []:
            self.addBin("<x> > -999")
        
        #print("Processing Bins...")
        for iBin,myBin in enumerate(self.bin_manager.bins):
            #print(f"...{iBin+1} of {len(self.bin_manager.bins)}...")
            myBin.sideband = copy.deepcopy(self.sideband); # Clone sideband structure
            
            if myBin.sideband["var"]=="": # User did not make a sideband region
                pars,errs=self.calculate_asymmetry(myBin.include_indecies)
                myBin.fitpars = pars
                myBin.fiterrs = errs
            else:
                self.calculate_asymmetries_for_each_sideband(myBin)
                
            
        # Print results
        if print_out:
            self.print()
        
        # Plot results
        if draw_asym:
            if self.bin_manager.bins[0].sideband["var"]!="":
                self.plot_asym_3()
            else:
                self.plot_asym()
            
        # Save to a yaml file
        if print_to_yaml:
            self.to_yaml(yaml_name)
            

class InjectionProject():
    def __init__(self , project_title = "", project_loc = "", n_trials = 10, n_cpus=4, injector_program_file = "",pion_pair="", infile = ""):
        self.project_title = project_title
        self.project_loc   = project_loc
        self.pion_pair     = pion_pair
        self.infile        = infile
        self.n_trials      = n_trials
        self.n_cpus        = n_cpus
        self.injector      = 0
        self.trial         = 0
        self.n_mods        = 0
        self.make_project_dir()
        self.seeds = np.random.randint(0,1000000,n_trials)
        self.yaml_files = []
        self.injector_program_file = injector_program_file
        self.load_injector_init_commands()
        
    def make_project_dir(self):
        print("="*100)
        if not os.path.exists(self.project_loc):
            os.makedirs(self.project_loc)
            print(f"Created directory `{self.project_loc}/`")
        else:
            print("Project Name",self.project_loc,"already taken...not recreating directory...")
        print("="*100)
    def load_injector_init_commands(self):

        # Add a prelude to the injector_init_commands to load the correct MC file
        injector_init_commands = f"self.injector = Injector('{self.pion_pair}','{self.infile}','{self.project_loc}')" + "\n" + self.injector_program_file.read()

        print("Injector Initialization Command:","\n","="*50,"\n",injector_init_commands,"\n\n")
        self.injector_init_commands = injector_init_commands
        self.n_mods        = self.count_mods()
    
    def create_injector(self):
        exec(self.injector_init_commands)

    def plot_modulation_result(self,mod, xlabel="",ylabel=""):
        
        yaml_files = self.yaml_files
        
        if self.has_sideband():
            data_path = "signal/y"
        else:
            data_path = "all/y"
        
        # Split path into components
        data_path_components = data_path.split('/')
        data_category, data_attribute = data_path_components

        plt.figure(figsize=(6,6),dpi=100)
        all_y_values = []
        true_y_values = None
        plt.axhline(0,color="black",linestyle="dashed")
        # Iterate over each YAML file
        for iyaml,yaml_file in tqdm(enumerate(yaml_files)):
            with open(yaml_file, 'r') as file:
                data = yaml.safe_load(file)

            mod_data = data[f"mod_{mod}"]
            x = data["x"]

            y_values = []
            y_errors = []
            # Iterate over each bin
            for bin_key in natsorted(mod_data.keys()):
                bin_data = mod_data[bin_key]
                if data_category in bin_data and data_attribute in bin_data[data_category]:
                    y_values.append(bin_data[data_category][data_attribute])
                    y_errors.append(bin_data[data_category][data_attribute+"err"])
            all_y_values.append(y_values)

        
            true_signal_y_values = []
            true_sigbkg_y_values = []
            true_bkg_y_values = []
            if "all" in data_path:
                true_paths = ["all/true_signal_y" , "all/true_bkg_y", "all/true_signal_y"]
            else:
                true_paths = ["signal/truey", "bkg/truey","signal+bkg/truey"]

            for true_y_values,true_path in zip([true_signal_y_values, true_bkg_y_values, true_sigbkg_y_values],true_paths):
                true_path_components = true_path.split('/')
                if len(true_path_components) != 2:
                    raise ValueError("true_path should be in format 'category/attribute'")
                true_category, true_attribute = true_path_components


                for bin_key in natsorted(mod_data.keys()):
                    bin_data = mod_data[bin_key]
                    if true_category in bin_data and true_attribute in bin_data[true_category]:
                        true_y_values.append(bin_data[true_category][true_attribute])

            if iyaml==0:
                plt.plot(x, y_values, color='grey', alpha=0.3)
                # An extra plot with no data, for a darker legend.
                plt.plot([], [], color='grey', alpha=0.3,label='Trials')
            elif iyaml==len(yaml_files)-1:
                plt.errorbar(x, y_values, yerr=y_errors,fmt='ko', capsize=3,label=f"Trial #{iyaml+1}")
            else:
                plt.plot(x, y_values, color='grey', alpha=0.3)
        # Convert to numpy array for easy calculations
        all_y_values = np.array(all_y_values)
        mean_y = np.mean(all_y_values, axis=0)
        std_y = np.std(all_y_values, axis=0)

        # Add error bands for 1σ and 3σ
        plt.fill_between(x, mean_y - std_y, mean_y + std_y, color='blue', alpha=0.3, label='1σ')
        plt.fill_between(x, mean_y - 3*std_y, mean_y + 3*std_y, color='red', alpha=0.2, label='3σ')

        # Plot the true y values if they were found
        plt.plot(x, true_signal_y_values, "r*", markersize=10,label='Injected (signal)')
        plt.plot(x, true_bkg_y_values, "r*", markersize=10,markerfacecolor="white",label='Injected (bkg)')
        plt.xlabel(data["xname"] if not xlabel else xlabel,fontsize=15)
        plt.ylabel(f"$A_{{LU}}^{{{mod_data['modname']}}}$",fontsize=15)
        ymax = np.amax(np.abs(plt.ylim()))
        plt.ylim(-ymax,ymax)

        plt.legend(fontsize=10)
        plt.savefig(self.project_loc+"/modulation_{}.png".format(mod))
        plt.close()
    
    
    def has_sideband(self, data=""):
        if data == "":
            with open(self.yaml_files[0], 'r') as file:
                data = yaml.safe_load(file)
            
        if "all" in data:
            return True
        for value in data.values():
            if isinstance(value, dict):
                if self.has_sideband(value):
                    return True
        return False
    
    def count_mods(self):
        with open(self.injector_program_file.name, 'r') as file:
            content = file.read()

        # Look for load_weight_funcs calls
        function_calls = re.findall(r'self\.injector\.load_weight_funcs\(\s*\[\s*(.*?)\s*\]\s*\)', content, re.DOTALL)

        # Initialize sublist counter
        sublist_count = 0

        # For each function call found
        for call in function_calls:
            # Count the sublists in the call
            sublists = re.findall(r'\[[^\[\]]*\]', call)
            sublist_count += len(sublists)

        return sublist_count


    def make_plots(self):
        
        self.yaml_files = [f"{self.project_loc}/trial_{i}.yaml" for i in range(self.n_trials)]
        
        # Save plots of the results    
        for i in range(self.n_mods):
            self.plot_modulation_result(mod = i)
            
        
            