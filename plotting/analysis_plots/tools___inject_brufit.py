import uproot
import uproot3 as uproot3
import numpy as np
import shutil
from array import array
import re
import ROOT
import matplotlib.pyplot as plt
import yaml
import os
import copy
import glob
from concurrent.futures import ProcessPoolExecutor
ROOT.gErrorIgnoreLevel = ROOT.kFatal
from natsort import natsorted # pip install natsort
from tqdm import tqdm
import array
from tools___io import *
from tools___etc import *


def get_azi_modulations(weight_funcs):
    
    """
    Function that returns a list of azimuthal modulation string for brufit
    """    
    
    _polarization = 1
    
    str_vec = []
    
    for i,wf in enumerate(weight_funcs):
        modulation = wf[2]
        modulation=modulation.replace(":phi","@phi").replace(":","[]")
        str_val = f"mod{i}={_polarization}*@hel[]*{modulation}"
        str_vec.append(str_val)
        
    return str_vec


class Bru:
    """
    BRUFIT class handler: Separated into two pieces
        One dedicated to handling the Pi0 sPlot
        One dedicated to handling the sWeighted azi fits
    """    
    def __init__(self, outdir, IDbranchname="fggID"):
        
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        
        self.outdir = outdir
        
        self.splot = ROOT.sPlot()
        self.splot.SetUp().SetOutDir(outdir+"/outSplot/")
        self.splot.SetUp().SetIDBranchName(IDbranchname)
        
        self.fm    = ROOT.FitManager()
        self.fm.SetUp().SetOutDir(outdir+"/outObsBins/")
        self.fm.SetUp().SetIDBranchName(IDbranchname)
        
        self.variables = []
        self.formulas  = []
        
    def load_splot_variable(self,v,vmin,vmax):
        
        self.splot.SetUp().LoadVariable(f"{v}[{vmin},{vmax}]")
        
    
    def load_splot_signal_pdf(self,s):
        
        self.splot.SetUp().FactoryPDF(s)
        self.splot.SetUp().LoadSpeciesPDF("Signal")    
    
    
    def load_splot_bkg_pdf(self,s):
        
        self.splot.SetUp().FactoryPDF(s)
        self.splot.SetUp().LoadSpeciesPDF("BG",1)
        
        
    def load_splot_data(self,rootfile,ttree):
        
        self.splot.LoadData(ttree,rootfile)
        
    
    def load_fm_bins_linspace(self,v,nbins,vmin,vmax):
        
        self.fm.Bins().LoadBinVar(v,nbins,vmin,vmax)
        
        
    def load_fm_bins(self,v,nbins,bins):
        
        self.fm.Bins().LoadBinVar(v,nbins,bins)
        
    
    def load_fm_variable(self,v,vmin,vmax):
        
        self.fm.SetUp().LoadVariable(f"{v}[{vmin},{vmax}]")
        self.variables.append(v)
    
    def load_fm_formula(self,formula):
        
        self.fm.SetUp().LoadFormula(formula)
        self.formulas.append(formula.split("=")[0])
        
        
    def load_fm_data(self,rootfile,ttree):
        
        self.fm.LoadData(ttree,rootfile)
        
        
    def load_fm_weights(self,weightfile):

        self.fm.Data().LoadWeights("Signal",weightfile)
        
        
    def run_splot(self):
        
        ROOT.Here.Go(self.splot)
        #self.splot.DrawWeighted("M2","Signal")
        #self.splot.DeleteWeightedTree()
    
    def run_fm(self):
        
        rhs = ":".join([f'b{i}[0,-1,1];{self.formulas[i]}' for i in range(len(self.formulas))])
        self.fm.SetUp().FactoryPDF(f"RooComponentsPDF::AziFit(1,{{{','.join(self.variables)}}},={rhs})")
        self.fm.SetUp().LoadSpeciesPDF("AziFit",1)
        
        
        ROOT.Here.Go(self.fm)
        
        
class BinManager:
    """
    BIN MANAGER class: Stores individual bin objects and carries out grouped functions
    """ 
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
    """
    BIN class: Individually stores bin/resulting fit info
    """
    def __init__(self,bin_func,include_indecies):
        self.bin_func = bin_func
        self.neat_bin_func = self.parse_bin_func(bin_func)
        self.include_indecies = include_indecies
        
        # The below two lists are populated when sidebands are not used
        self.fitpars  = []
        self.fiterrs  = [] 
    
    @staticmethod
    def parse_bin_func(bf):
        pattern = r':(\w+):'
        return re.sub(pattern, r'\1', bf)
    
    

class Injector:
    ####################################################################################
    # Intializer + Initializer Functions
    ####################################################################################
    def __init__(self, pion_pair, filename):
        self.filename = filename
        self.u = uproot.open(self.filename)
        self.u = self.u["dihadron_cuts"]
        self.N = self.u.num_entries
        self.M2 = self.u["M2"].array(library="np")
        self.pion_pair = pion_pair
        self.set_pion_pids()
        self.load_azi_vars()
        self.weight_funcs = [ ["0.2","0.0","sin(:phi_R0:)"] ]
        self.bPol = 1.0
        self.nbins = 30 # For 2d fits

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
    # ex: weight_funcs = [ ["0.05*<Mh>" , "0.10*<Mh>", "sin(<phi_h>-<phi_R0>)"],
    #                      ["-0.1*<z>" , "0.0*<z>"] , "sin(<phi_h>)"] 
    #                    ]
    def load_weight_funcs(self, weight_funcs):
        self.weight_funcs = weight_funcs
    
    # Set the beam polarization [0,1]
    def set_beam_polarization(self, bPol):
        self.bPol = bPol

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
            
    def set_outdir(self,outdir):
        # Make the directory in volatile so that we can save our trees to it
        trueoutdir = "/work/clas12/users/gmat/clas12/clas12_dihadrons/projects/pipi0_paper_RGA_only/volatile/asym/inject_bru/"+outdir
        if not os.path.exists(trueoutdir): 
            os.makedirs(trueoutdir)
        self.outdir = os.getcwd()+"/output_plots/pipi0_paper_RGA_only/inject_bru/"+outdir
        try:
            os.symlink(trueoutdir,self.outdir)
        except:
            pass
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
    
    ####################################################################################
    # Save to yaml file
    ####################################################################################  
    def to_yaml(self,yaml_name = "out.yaml"):
        
        
        data = self.bin_manager.get_bin_avgs()[0] # Get bin centers (x-values)
        x = data["x"]
        outyaml = {"x":[float(X) for X in x],
                   "xname": self.bin_manager.bin_names[0]}
        for iwf,wf in enumerate(self.weight_funcs):
            outyaml[f"mod_{iwf}"] = {"modname":self.parse_func_latex(wf[2])}
            asym_sig_values, asym_bkg_values = self.get_asym_values(iwf)
            
            for iBin,myBin in enumerate(self.bin_manager.bins):
                u = uproot.open(self.outdir+f"/bin_{iBin}/outObsBins/ResultsHSMinuit2.root")
                u = u["ResultTree"]
                
                outyaml[f"mod_{iwf}"][f"bin_{iBin}"] = {"all":{}}
                outyaml[f"mod_{iwf}"][f"bin_{iBin}"]["all"]["y"] = float(u[f"b{iwf}"].array(library="np")[0])
                outyaml[f"mod_{iwf}"][f"bin_{iBin}"]["all"]["yerr"] = float(u[f"b{iwf}_err"].array(library="np")[0])
                outyaml[f"mod_{iwf}"][f"bin_{iBin}"]["all"]["true_signal_y"] = float(myBin.get_bin_avg_of_branch(asym_sig_values))
                outyaml[f"mod_{iwf}"][f"bin_{iBin}"]["all"]["true_bkg_y"] = float(myBin.get_bin_avg_of_branch(asym_bkg_values))
        
        with open(yaml_name, 'w') as file:
            yaml.dump(outyaml, file, default_flow_style=False)
            
    
    ####################################################################################
    # Execution command
    ####################################################################################            
    # Calculate the asymmetry histograms
    def run(self,yaml_file=""):
        # Produce helicities by weighting procedure
        self.get_helicities() 
        
        # Save binned TTrees
        _M2    = array.array('d',[0.0])
        _phi_h = array.array('d',[0.0])
        _phi_R = array.array('d',[0.0])
        _hel   = array.array('i',[0])
        _fggID   = array.array('i',[0])
        
        
        tfilename = self.outdir+"/tmp.root"
        print(tfilename)
        self.tfilename = tfilename

        with uproot3.recreate(tfilename) as fOut:

            for iBin, myBin in tqdm(enumerate(self.bin_manager.bins)):
                # Create a mask array using vectorized operations
                mask = np.array(myBin.include_indecies) != 0
                # Filter data using the mask
                PhiH_bin = np.array(self.PhiH)[mask]
                PhiR_bin = np.array(self.PhiR0)[mask]
                hel_bin = np.array(self.hel)[mask]
                M2_bin = np.array(self.M2)[mask]
                # Create fggID for the bin
                fggID_bin = np.arange(len(PhiH_bin))

                # Create tree dictionary
                data = {
                    "phi_h": PhiH_bin,
                    "phi_R0": PhiR_bin,
                    "hel": hel_bin.astype(np.int32),
                    "M2": M2_bin,
                    "fggID": fggID_bin
                }

                # Define dtypes for each branch
                dtypes = {
                    "phi_h": np.float64,
                    "phi_R0": np.float64,
                    "hel": np.int32,
                    "M2": np.float64,
                    "fggID": np.int32
                }

                # Create and fill the tree
                fOut[f"bin_{iBin}"] = uproot3.newtree(dtypes)
                fOut[f"bin_{iBin}"].extend(data)
                

        for iBin,myBin in tqdm(enumerate(self.bin_manager.bins)):
            self.run_bru_on_bin(iBin,tfilename)
            
        if yaml_file!="":
            self.to_yaml(yaml_file)
        
    def run_bru_on_bin(self,iBin,tfilename):

        myBin = self.bin_manager.bins[iBin]
        u = uproot.open(tfilename)
        u = u[f"bin_{iBin}"]
        M2max = np.amin([0.4,np.amax(u["M2"].array(library="np"))])
        bru = Bru(outdir = self.outdir+f"/bin_{iBin}")
        bru.load_splot_variable("M2",0.03,M2max)
        bru.load_splot_signal_pdf("Gaussian::Signal( M2, mean[0.131,0.129,0.140], sigma[0.01,0.001,0.02])")
        bru.load_splot_bkg_pdf("Chebychev::BG(M2,{a0[0,-1,1],a1[0,-1,1],a2[0,-1,1]})")
        bru.load_splot_data(tfilename,f"bin_{iBin}")
        bru.run_splot()

        bru.load_fm_variable("phi_h",-3.1415,3.1415)
        bru.load_fm_variable("phi_R0",-3.1415,3.1415)
        bru.load_fm_variable("hel",-2,2)
        for formula in get_azi_modulations(self.weight_funcs):
            bru.load_fm_formula(formula)
        #bru.load_fm_bins_linspace("Mh",4,0.3,1.3)
        bru.load_fm_data(tfilename,f"bin_{iBin}")
        bru.load_fm_weights(self.outdir+f"/bin_{iBin}/outSplot/Weights.root")
        bru.run_fm()
        
        

class InjectionBruProject():
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
        injector_init_commands = f"self.injector = Injector('{self.pion_pair}','{self.infile}')" + "\n"
        injector_init_commands += self.injector_program_file.read()

        print("Injector Initialization Command:","\n","="*50,"\n",injector_init_commands,"\n\n")
        self.injector_init_commands = injector_init_commands
        self.n_mods        = self.count_mods()
    
    def create_injector(self,trial_num):
        self.injector_init_commands += "\n" + f"self.injector.set_outdir('{self.project_title}/trial_{trial_num}')"
        exec(self.injector_init_commands)
        
    def plot_modulation_result(self,mod, xlabel="",ylabel=""):
        
        yaml_files = self.yaml_files
        
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
            true_paths = ["all/true_signal_y" , "all/true_bkg_y", "all/true_signal_y"]
            
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