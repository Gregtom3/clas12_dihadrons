import ROOT
from binning import *
import copy
import matplotlib.pyplot as plt
from collections import defaultdict
from array import array
from utils import *
import math 

def create_1d_boundary_tgraph(boundary_exprs, n_points=1000, repl=[], yrange=(0.001,1000000)):
    # Split the region into individual segments
    boundary_exprs = boundary_exprs.split("&")
    # Get the min and max from the boundary exprs
    xmin = extract_number_with_decimal(boundary_exprs[0])
    xmax = extract_number_with_decimal(boundary_exprs[1])
    # Calculate an array for x,y
    xmin = list(xmin*np.ones(n_points))
    xmax = list(xmax*np.ones(n_points))
    y = np.linspace(*yrange,n_points)
    # Create the left and right TGraph boundary
    tgraphs = []
    tgraphs.append(ROOT.TGraph(n_points,array('d',xmin),array('d',y)))
    tgraphs.append(ROOT.TGraph(n_points,array('d',xmax),array('d',y)))
    xtitle = retitle(repl[0])
    tgraphs[0].SetTitle(f";{xtitle};")
    tgraphs[1].SetTitle(f";{xtitle};")
    
    return copy.deepcopy(tgraphs)
    
def create_2d_boundary_tgraph(boundary_exprs, n_points=1000,repl=["x","Q2"]):
    # Split the region into individual segments
    boundary_exprs = boundary_exprs.split("&")
    # Define the range of x and y values
    if(repl==["x","Q2"]):
        x_vals = np.linspace(0, 1, 1000)
        y_vals = np.linspace(0, 20, 1000)
    elif(repl==["z","pTtot"]):
        x_vals = np.linspace(-0.1, 1.1, 1000)
        y_vals = np.linspace(-0.1, 3, 1000)
    else:
        x_vals = np.linspace(0, 1, 1000)
        y_vals = np.linspace(0, 20, 1000)
    # Create a meshgrid from x and y values
    x, y = np.meshgrid(x_vals, y_vals)

    # Evaluate the inequalities for each point in the meshgrid
    ineq = np.ones_like(x, dtype=bool)
    for inequality in boundary_exprs:
        inequality=inequality.replace(repl[1],"y").replace(repl[0],"x")

        ineq = np.logical_and(ineq, eval(inequality))

    # Define the levels of the contour plot
    levels = [0.5, 1]

    # Create a contour plot of the inequalities
    contour_data= plt.contour(x, y, ineq, levels=levels, colors='k')
    # Hide the plot
    plt.close()
    xy_points = np.array(contour_data.collections[0].get_paths()[0].vertices)
    x=xy_points[:,0]
    y=xy_points[:,1]
    
    # Create the TGraph object
    graph = ROOT.TGraph(len(x), array('d', x), array('d', y))

    # Set the title and axis labels
    graph.SetTitle("")
    xtitle = retitle(repl[0])
    ytitle = retitle(repl[1])
    graph.GetXaxis().SetTitle(xtitle)
    graph.GetYaxis().SetTitle(ytitle)
    graph.SetLineColor(1)
    graph.SetLineWidth(2)
    
    return graph

def get_tgraph_boundaries(input_boundaries,repl=["x","Q2"]):
    btgraphs=[]
    input_boundaries=remove_duplicates(input_boundaries)
    
    for b in input_boundaries:
        btgraph = (create_2d_boundary_tgraph(b,1000,repl) if len(repl)==2 else create_1d_boundary_tgraph(b,repl=repl))
        if(type(btgraph)==list): # 1d_boundary returns a list of TGraphs
            for bt in btgraph:
                btgraphs.append(bt)
        else:
            btgraphs.append(btgraph)
    return copy.deepcopy(btgraphs)

def get_graph_bounds(tgraphs):
    xmin = ymin = float('inf')
    xmax = ymax = float('-inf')

    for i in range(len(tgraphs)):
        x_min_graph = tgraphs[i].GetXaxis().GetXmin()
        x_max_graph = tgraphs[i].GetXaxis().GetXmax()
        y_min_graph = tgraphs[i].GetYaxis().GetXmin()
        y_max_graph = tgraphs[i].GetYaxis().GetXmax()

        if x_min_graph < xmin:
            xmin = x_min_graph
        if x_max_graph > xmax:
            xmax = x_max_graph
        if y_min_graph < ymin:
            ymin = y_min_graph
        if y_max_graph > ymax:
            ymax = y_max_graph

    return xmin, xmax, ymin, ymax

def set_bounds(tgraphs):
    
    xmin, xmax, ymin, ymax = get_graph_bounds(tgraphs)
    tgraphs[0].GetXaxis().SetLimits(xmin, xmax)
    tgraphs[0].GetYaxis().SetRangeUser(ymin, ymax)
    

# def plot_map(bin_manager, hists_pars=None, data_tree=None, mc_tree=None, uID=None, plot_variable=None, extra_cut=None):
#     # Extract binning information
#     rect_bins = bin_manager.rect_bins
#     rect_names = bin_manager.rect_names
#     custom_bins = bin_manager.custom_bins
#     custom_names = bin_manager.custom_names

#     # Format bins and names for plotting
#     rect_bins = reformat_rect_bins(rect_bins, rect_names)
#     rect_bins = split_into_fours(rect_bins)
#     rect_names = split_into_twos(rect_names)
#     custom_bins = [custom_bin.replace("[","").replace("]","") for custom_bin in custom_bins]

#     # Get the count of canvas pads
#     custom_pad_count = 1 if bin_manager.has_custom_bins else 0
#     rect_pad_count = len(rect_bins) if bin_manager.has_rect_bins else 0
#     hist_pad_count = 1 if plot_variable else 0
    
#     total_pads = custom_pad_count + rect_pad_count + hist_pad_count
#     # Make sure we have enough user defined histograms
#     if(len(hists_pars)!=total_pads):
#         print("ERROR: User must provide a list of histograms to be plotted into...Aborting...")
#     # Create the hists
#     hists=[]
#     for i in range(total_pads):
#         hists.append(ROOT.TH1F(*hists_pars[i]) if len(hists_pars[i])==5 else ROOT.TH2F(*hists_pars[i]))
#         hists[i].SetName(f"hist_{i}")
#     if mc_tree!=None:
#         hists.append(ROOT.TH1F(*hists_pars[-1]) if len(hists_pars[-1])==5 else ROOT.TH2F(*hists_pars[-1]))
#         hists[i+1].SetName(f"hist_{i+1}")
        
#     # Get the unique custom and rect bin index from uID
#     if uID is not None:
#         rect_index, custom_index = bin_manager.get_bin_ids_from_unique_id(uID)

#     ROOT.gStyle.SetOptStat(0)
#     # Set histogram line width
#     ROOT.gStyle.SetHistLineWidth(2)

#     # Set left and right margins
#     ROOT.gStyle.SetPadLeftMargin(0.15)
#     ROOT.gStyle.SetPadRightMargin(0.15)

#     # Set axes title size
#     ROOT.gStyle.SetTitleSize(0.05, "X")  # X-axis title size
#     ROOT.gStyle.SetTitleSize(0.05, "Y")  # Y-axis title size
    
#     # Prepare canvas
#     canvas_height = 600
#     canvas_width = 600
#     canvas = ROOT.TCanvas("c","c", canvas_width * total_pads, canvas_height)
#     canvas.Divide(total_pads,1)
#     hist_num = 0
#     # Draw custom binning pad
#     if bin_manager.has_custom_bins:
#         canvas.cd(1)
#         custom_tgraphs = get_tgraph_boundaries(custom_bins, ["rec_x", "rec_Q2"])

#         if data_tree is not None:
#             data_tree.Draw(f"{custom_names[1]}:{custom_names[0]}>>hist_{hist_num}", "rec_passDihadron", "colz")
#             hist_num+=1
#             ROOT.gPad.SetLogz(1)

#         # Draw boundary for each custom bin
#         for i, graph in enumerate(custom_tgraphs):
#             draw_option = "L same" if i > 0 else ("A" if data_tree is None else "") + "L"
#             graph.Draw(draw_option)

#         # Highlight selected bin
#         if uID is not None:
#             custom_tgraphs[custom_index].SetLineColor(2)
#             custom_tgraphs[custom_index].SetLineWidth(3)
#             custom_tgraphs[custom_index].Draw("L same")

#         set_bounds(custom_tgraphs)

#     # Draw rectangular binning pad
#     if bin_manager.has_rect_bins:
#         rect_tgraphs = [get_tgraph_boundaries(rect_bin, rect_name) for rect_bin, rect_name in zip(rect_bins, rect_names)]

#         for pad_num in range(rect_pad_count):
#             canvas.cd(pad_num + 1 + custom_pad_count)
#             rect_tgraphs_current_pad = rect_tgraphs[pad_num]

#             if data_tree is not None:
#                 # Apply a cut on the drawn distribution based on the panels before
#                 full_tcut = ["rec_passDihadron"]
#                 if(bin_manager.has_custom_bins):
#                     full_tcut.append(custom_bins[custom_index])
#                 for rect_bin_num in range(pad_num-1):
#                     full_tcut.append(rect_bins[rect_bin_num][rect_index])
#                 full_tcut = "&".join(full_tcut)
#                 full_tcut = full_tcut.replace("&","&&")
#                 if len(rect_names[pad_num]) == 2:
#                     data_tree.Draw(f"{rect_names[pad_num][1]}:{rect_names[pad_num][0]}>>hist_{hist_num}", full_tcut, "colz")
#                 else:
#                     data_tree.Draw(f"{rect_names[pad_num][0]}>>hist_{hist_num}", full_tcut, "hist")
#                 hist_num+=1
                
#             # Draw boundary for each rectangular bin
#             for i, graph in enumerate(rect_tgraphs_current_pad):
#                 draw_option = "L same" if i > 0 else ("A" if data_tree is None else "") + "L"
#                 graph.Draw(draw_option)

#             # Highlight selected bin
#             if uID is not None:
#                 rect_tgraphs_current_pad[rect_index].SetLineColor(2)
#                 rect_tgraphs_current_pad[rect_index].SetLineWidth(3)
#                 rect_tgraphs_current_pad[rect_index].Draw("L same")

#             set_bounds(rect_tgraphs_current_pad)

#     # Draw histogram
#     if plot_variable!=None:
#         # cd into the last pad
#         canvas.cd(total_pads)
#         # Get full cut for the histogram
#         bin_tcut = bin_manager.get_tcut_from_unique_id(uID)
#         full_tcut = bin_tcut if extra_cut==None else "&&".join([bin_tcut,extra_cut])
#         # Draw into the user provided histogram
#         data_tree.Draw(f"{plot_variable}>>hist_{hist_num}",full_tcut,"hist")
#         hists[hist_num].SetLineColor(ROOT.kBlack);
#         hist_num+=1;
#         if mc_tree!=None:
#             hists[hist_num].SetLineColor(ROOT.kBlue);
#             mc_tree.Draw(f"{plot_variable}>>hist_{hist_num}",full_tcut,"hist same")
#             hists[hist_num].Scale(data_tree.GetEntries("rec_passDihadron")/mc_tree.GetEntries("rec_passDihadron"))
            
#             # Set ymax
#             ymax = np.amax([hists[hist_num-1].GetMaximum(),hists[hist_num].GetMaximum()])
#             hists[hist_num-1].GetYaxis().SetRangeUser(0,ymax*1.205)
#             # Create TLegend
#             legend = ROOT.TLegend(0.15, 0.8, 0.85, 0.9)
#             legend.SetNColumns(2)
#             legend.AddEntry(hists[hist_num - 1], "Data", "l")
#             legend.AddEntry(hists[hist_num], "Monte Carlo", "l")
#             legend.Draw()
            
#             # Set the Intersection Over Union as title
#             iou = calculate_iou(hists[hist_num-1],hists[hist_num])
#             hists[hist_num-1].SetTitle(f"Intersection over Union = {iou:.4f}")
#     return copy.deepcopy(canvas)


def plot_map(bin_manager, histograms, data_tree, mc_tree, uID=None):
    # Extract binning information
    rect_bins = bin_manager.rect_bins
    rect_names = bin_manager.rect_names
    custom_bins = bin_manager.custom_bins
    custom_names = bin_manager.custom_names

    # Format bins and names for plotting
    rect_bins = reformat_rect_bins(rect_bins, rect_names)
    rect_bins = split_into_fours(rect_bins)
    rect_names = split_into_twos(rect_names)
    custom_bins = [custom_bin.replace("[","").replace("]","") for custom_bin in custom_bins]

    # Get the count of canvas pads
    custom_pad_count = 1 if bin_manager.has_custom_bins else 0
    rect_pad_count = len(rect_bins) if bin_manager.has_rect_bins else 0
    hist_pad_count = 1 
    
    total_pads = custom_pad_count + rect_pad_count + hist_pad_count
    
    rect_index, custom_index = bin_manager.get_bin_ids_from_unique_id(uID)

    ROOT.gStyle.SetOptStat(0)
    # Set histogram line width
    ROOT.gStyle.SetHistLineWidth(2)

    # Set left and right margins
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadRightMargin(0.15)

    # Set axes title size
    ROOT.gStyle.SetTitleSize(0.05, "X")  # X-axis title size
    ROOT.gStyle.SetTitleSize(0.05, "Y")  # Y-axis title size
    
    # Prepare canvas
    canvas_height = 600
    canvas_width = 600
    canvas = ROOT.TCanvas("c","c", canvas_width * total_pads, canvas_height)
    canvas.Divide(total_pads,1)
    hist_num = 0
    # Draw custom binning pad
    if bin_manager.has_custom_bins:
        canvas.cd(1)
        custom_tgraphs = get_tgraph_boundaries(custom_bins, ["rec_x", "rec_Q2"])
        
        histograms[0][uID].Draw("colz")
        hist_num+=1
        ROOT.gPad.SetLogz(1)

        # Draw boundary for each custom bin
        for i, graph in enumerate(custom_tgraphs):
            draw_option = "L same" if i > 0 else ("A" if data_tree is None else "") + "L"
            graph.Draw(draw_option)

        # Highlight selected bin
        if uID is not None:
            custom_tgraphs[custom_index].SetLineColor(2)
            custom_tgraphs[custom_index].SetLineWidth(3)
            custom_tgraphs[custom_index].Draw("L same")

        set_bounds(custom_tgraphs)

    # Draw rectangular binning pad
    if bin_manager.has_rect_bins:
        rect_tgraphs = [get_tgraph_boundaries(rect_bin, rect_name) for rect_bin, rect_name in zip(rect_bins, rect_names)]

        for pad_num in range(rect_pad_count):
            canvas.cd(pad_num + 1 + custom_pad_count)
            rect_tgraphs_current_pad = rect_tgraphs[pad_num]
            histograms[pad_num+custom_pad_count][uID].Draw("hist colz")
            hist_num+=1
            
            # Draw boundary for each rectangular bin
            for i, graph in enumerate(rect_tgraphs_current_pad):
                draw_option = "L same" if i > 0 else ("A" if data_tree is None else "") + "L"
                graph.Draw(draw_option)

            # Highlight selected bin
            if uID is not None:
                rect_tgraphs_current_pad[rect_index].SetLineColor(2)
                rect_tgraphs_current_pad[rect_index].SetLineWidth(3)
                rect_tgraphs_current_pad[rect_index].Draw("L same")

            set_bounds(rect_tgraphs_current_pad)

    # cd into the last pad
    canvas.cd(total_pads)

    histograms[-2][uID].SetLineColor(ROOT.kBlack)
    histograms[-1][uID].SetLineColor(ROOT.kBlue)

    histograms[-2][uID].Draw("hist")
    histograms[-1][uID].Draw("hist same")
    histograms[-1][uID].Scale(data_tree.GetEntries("rec_passDihadron")/mc_tree.GetEntries("rec_passDihadron"))

    # Set ymax
    ymax = np.amax([histograms[-1][uID].GetMaximum(),histograms[-2][uID].GetMaximum()])
    histograms[-2][uID].GetYaxis().SetRangeUser(0,ymax*1.205)
    # Create TLegend
    legend = ROOT.TLegend(0.15, 0.8, 0.85, 0.9)
    legend.SetNColumns(2)
    legend.AddEntry(histograms[-2][uID].GetValue(), "Data", "l")
    legend.AddEntry(histograms[-1][uID].GetValue(), "Monte Carlo", "l")
    legend.Draw()

    # Set the Intersection Over Union as title
    iou = calculate_iou(histograms[-2][uID],histograms[-1][uID])
    histograms[-2][uID].SetTitle(f"Intersection over Union = {iou:.4f}")
    return copy.deepcopy(canvas)


def make_multiple_bin_histos(bin_manager, data_tree, mc_tree, hists_pars=None, plot_variable = None, extra_cut = "rec_passDihadron"):
    
    df = ROOT.RDataFrame(data_tree)
    tdf = ROOT.RDataFrame(mc_tree)
    
    nbins = bin_manager.total_bins
    # Extract binning information
    rect_bins = bin_manager.rect_bins
    rect_names = bin_manager.rect_names
    custom_bins = bin_manager.custom_bins
    custom_names = bin_manager.custom_names

    # Format bins and names for plotting
    rect_bins = reformat_rect_bins(rect_bins, rect_names)
    rect_bins = split_into_fours(rect_bins)
    rect_names = split_into_twos(rect_names)
    custom_bins = [custom_bin.replace("[","").replace("]","") for custom_bin in custom_bins]

    # Get the count of canvas pads
    custom_pad_count = 1 if bin_manager.has_custom_bins else 0
    rect_pad_count = len(rect_bins) if bin_manager.has_rect_bins else 0
    total_pads = custom_pad_count+rect_pad_count    

    histograms = [[] for _ in range(custom_pad_count + rect_pad_count + 2)]
    
    for pad_num in range(total_pads):
        
        # Assuming df is your RDataFrame, the bins for each uID are stored in some way
        for uID in range(1,nbins):
            # Apply filters based on uID
            rect_index, custom_index = bin_manager.get_bin_ids_from_unique_id(uID)
            # Apply a cut on the drawn distribution based on the panels before
            full_tcut = ["rec_passDihadron"]
            if(bin_manager.has_custom_bins and pad_num>=1):
                full_tcut.append(custom_bins[custom_index])
            for rect_bin_num in range(bin_manager.has_custom_bins,pad_num):
                full_tcut.append(rect_bins[rect_bin_num][rect_index])
            full_tcut = "&".join(full_tcut)
            full_tcut = full_tcut.replace("&","&&")
            if bin_manager.has_custom_bins:
                if pad_num == 0:
                    names = custom_names
                else:
                    names = rect_names[pad_num-1]
            else:
                names = rect_names[pad_num]
            
            if(len(hists_pars[pad_num])==4): # 1D Histo
                hist = df.Filter(full_tcut).Histo1D((f"hist_{uID}_{pad_num}",*hists_pars[pad_num]),names[0])
                histograms[pad_num].append(hist)
            else:
                
                hist = df.Filter(full_tcut).Histo2D((f"hist_{uID}_{pad_num}",*hists_pars[pad_num]),names[0],names[1])
                histograms[pad_num].append(hist)
    
            # Get full cut for the histogram
            bin_tcut = bin_manager.get_tcut_from_unique_id(uID)
            full_tcut = bin_tcut if extra_cut==None else "&&".join([bin_tcut,extra_cut])
            # Draw into the user provided histogram
            rechist = df.Filter(full_tcut).Histo1D((f"hist_{uID}_{total_pads}",*hists_pars[total_pads]),plot_variable)
            
            histograms[total_pads].append(rechist)
            mchist = tdf.Filter(full_tcut).Histo1D((f"hist_{uID}_{total_pads+1}",*hists_pars[total_pads]),plot_variable)
            histograms[total_pads+1].append(mchist)
            
    return histograms
    
    
    
    


class MultidimensionalHist:
    def __init__(self,bin_manager=None, histo=None):
        self.labels = []
        self.bins = defaultdict(lambda: defaultdict(list))
        if histo is not None and bin_manager is not None:
            self.load_hist(bin_manager,histo)
        
    def load_hist(self,bin_manager,histo):

        rect_names = bin_manager.rect_names
        custom_names = bin_manager.custom_names # List like ["rec_x","rec_Q2"]
        custom_names = "-".join([cn.replace("rec_","") for cn in custom_names]) # --> ["x-Q2"]

        bin_names = (rect_names if bin_manager.has_rect_bins else []) + ([custom_names] if bin_manager.has_custom_bins else [])
        self.labels = bin_names

        for uID in range(1,histo.GetNbinsX()+1):
            rect_id, custom_id = bin_manager.get_bin_ids_from_unique_id(uID)
            rect_ids = bin_manager.get_all_rect_bin_ids_from_rect_bin_id(rect_id)
            rect_bins = bin_manager.rect_bins
            rect_lefts = np.array([rect_bins[i][rect_id] for i,rect_id in enumerate(rect_ids)])
            rect_rights = np.array([rect_bins[i][rect_id+1] for i,rect_id in enumerate(rect_ids)])
            rect_centers = 0.5 * (rect_lefts+rect_rights)
            rect_centers = np.round(rect_centers,4)

            ids = (rect_ids if bin_manager.has_rect_bins else []) + ([custom_id] if bin_manager.has_custom_bins else [])
            centers = (list(rect_centers) if bin_manager.has_rect_bins else []) + ([int(custom_id)] if bin_manager.has_custom_bins else [])
            lefts = (list(rect_lefts) if bin_manager.has_rect_bins else []) + ([custom_id-0.5] if bin_manager.has_custom_bins else [])
            rights = (list(rect_rights) if bin_manager.has_rect_bins else []) + ([custom_id+0.5] if bin_manager.has_custom_bins else [])
            counts = histo.GetBinContent(uID)
            error = histo.GetBinError(uID)

            ids = tuple(ids)
            centers = np.array(centers)
            lefts = np.array(lefts)
            rights = np.array(rights)

            self.add_data(ids, lefts, centers, rights, counts, error)
            
    def add_data(self, indices, left, center, right, data, error):
        self.bins[indices]['center'].append(center)
        self.bins[indices]['data'].append(data)
        self.bins[indices]['error'].append(error)
        self.bins[indices]['left'].append(left)
        self.bins[indices]['right'].append(right)

    def plot(self, *args, plot_type="bar", make_title=True, log_scale=False):
        fixed_indices = [idx for idx, arg in enumerate(args) if type(arg) is int]
        xaxis_index = next((idx for idx, arg in enumerate(args) if arg == "<X-AXIS>"), None)
        labels_index = next((idx for idx, arg in enumerate(args) if arg == "<LABELS>"), None)

        plt.figure()
        plot_data = defaultdict(lambda: {'x': [], 'y': [], 'xerr': [], 'width': []})

        for indices, bin in self.bins.items():
            if all(indices[fixed_idx] == args[fixed_idx] for fixed_idx in fixed_indices):
                the_center = self.bins[indices]['center'][0][labels_index] if labels_index is not None else None
                label = f"{retitle(self.labels[labels_index])}={the_center}" if labels_index is not None else None
                plot_data[label]['x'].extend([center[xaxis_index] for center in bin['center']])
                plot_data[label]['y'].extend(bin['data'])
                plot_data[label]['xerr'].extend(bin['error'])
                plot_data[label]['width'].extend([right[xaxis_index]-left[xaxis_index] for right,left in zip(bin['right'],bin['left'])])
                
        title=""
        if make_title is True and fixed_indices:
            title_parts = []
            for fixed_idx in fixed_indices:
                label = retitle(self.labels[fixed_idx])
                converted_args = tuple(0 if isinstance(item, str) else item for item in tuple(args))
                center_value = self.bins[converted_args]["center"][0][fixed_idx]
                if("-" in label): # True for custom bins
                    title_parts.append(f"${label}$ bin {int(center_value)}")
                else:
                    title_parts.append(f"${label}$={center_value}")
            title = ", ".join(title_parts)

        draw_legend=True

        for label, data in plot_data.items():
            if label is None:
                draw_legend=False
            if plot_type == "bar":
                plt.bar(data['x'], data['y'], yerr=data['xerr'], width=data['width'], label=label,capsize=3,ec="black",alpha=0.6)
            elif plot_type == "line":
                plt.errorbar(data['x'], data['y'], yerr=data['xerr'], fmt='o-', label=label)

        xtitle = retitle(self.labels[xaxis_index]) if xaxis_index is not None else "Index" 
        plt.xlabel(f"${xtitle}$",fontsize=15)
        plt.ylabel('Counts',fontsize=15)
        if log_scale:
            plt.yscale('log')
        if title:
            plt.title(title)
        if draw_legend:
            ncol = math.ceil(len(plot_data) / 5)  # Ensure at most 5 rows
            plt.legend(ncol=ncol, loc="upper center")  # Adjust position to avoid overlapping

        plt.show()

    def print_usage(self):
        print("Usage for the plot method:")
        print("plot(*args, plot_type='bar', make_title=True)")
        print("args: Integers or '<X-AXIS>', '<LABELS>' representing dimensions. Use integers to fix a dimension to a specific bin.")
        print("plot_type: 'bar' or 'line'. Default is 'bar'.")
        print("make_title: Boolean. If True, a title will be generated from the fixed indices. Default is True.")
        print("log_scale: Boolean. If True, the y-axis will be set to a logarithmic scale. Default is False.")
        
        print("\nDimensions and their labels:")
        for i, label in enumerate(self.labels):
            print(f"Dimension {i}: {label}")

        print("\nNumber of bins in each dimension:")
        for i in range(len(self.labels)):
            unique_bins_in_dim = len(set(key[i] for key in self.bins.keys()))
            print(f"Dimension {i}: {unique_bins_in_dim} bins")
        