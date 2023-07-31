import os
import numpy as np
import uproot
import yaml
import re
import math
import ROOT
import copy
import concurrent.futures
from tqdm import tqdm

def print_yaml_structure(data, indent=0):
    for key, value in data.items():
        print(" " * indent + key + ":")
        if isinstance(value, dict):
            print_yaml_structure(value, indent+2)
        else:
            if(type(value)==list):
                value="<LIST>"
            print(" " * (indent+2) + str(value))
            
def extract_numbers(string):
    pattern = r"-?\d+(?:\.\d+)?"
    numbers = re.findall(pattern, string)
    return [float(number) for number in numbers]

def get_directories_with_values(directory_path, x_values):
    # get a list of directories in the given path
    directories = os.listdir(directory_path)
    
    # filter the list to only include directories that contain every element in x_values
    filtered_directories = []
    for d in directories:
        # Sometimes, brufit rounds the bin weirdly, so I do this to make sure we stil get it
        # Be careful though if the binning is really tiny!
        if (all(str(x) in d for x in x_values)) or (all(str(np.round(x-0.01,2)) in d for x in x_values)) or (all(str(np.round(x+0.01,2)) in d for x in x_values)):
            filtered_directories.append(d)
    
    return filtered_directories

def read_bruout(result_file):
    param_names = []
    param_values = []
    param_errors = []
    
    with uproot.open(result_file) as f:
        if 'ResultTree' in f:
            tree = f['ResultTree']

            # Loop over all branches in the tree
            for branch_name in tree.keys():

                # Skip the error branches
                if branch_name.endswith('_err'):
                    continue
                if(branch_name in ["NLL"]):
                    continue
                # Get the parameter name and its value
                param_name = branch_name
                param_value = tree[param_name].array()[0]

                # Get the corresponding error branch
                err_branch_name = f'{param_name}_err'
                if err_branch_name not in tree:
                    continue

                # Get the parameter error
                param_error = tree[err_branch_name].array()[0]

                # Append the data to the arrays
                if param_name not in param_names:
                    param_names.append(param_name)
                    param_values.append([float(param_value)])
                    param_errors.append([float(param_error)])
                else:
                    index = param_names.index(param_name)
                    param_values[index].append(float(param_value))
                    param_errors[index].append(float(param_error))
    return param_names, param_values, param_errors

def read_splot(subdir):
    relpath = subdir.split('/')[-1]
    
    x_values = extract_numbers(relpath)

    # Open the TTree called ResultTree
    result_file = os.path.join(subdir, 'ResultsHSMinuit2.root')
    if(not os.path.exists(result_file)):
        return None,None,None,None
    param_names, param_values, param_errors = read_bruout(result_file)
    return x_values , param_names, param_values, param_errors

def read_sideband(subdir=""):
    relpath = subdir.split('/')[-1]
    
    x_values = extract_numbers(relpath)
    
    def calculate_means(lst):
        result = []
        for i in range(0, len(lst), 2):
            mean = (lst[i] + lst[i+1]) / 2
            result.append(np.round(mean,2))
        return result
    
    x_values = calculate_means(x_values)
    
    # Open the TTree
    result_file = os.path.join(subdir, 'sideband.root')
    if(not os.path.exists(result_file)):
        return None,None,None,None
     # open the result file with PyROOT
    root_file = ROOT.TFile(result_file)
    # access the TVectorT<double> object
    purity_vector = root_file.Get("purity_4")
    purity = purity_vector[0]
    root_file.Close()
    
    # Now get the sigbg sideband data
    sigbg_file="/".join(subdir.split('/')[:-2])+"/outObsBins_sdbnd_sigbg"
    sigbg_file+="/"+get_directories_with_values(sigbg_file, x_values)[0]+"/ResultsHSMinuit2.root"
    if(not os.path.exists(sigbg_file)):
        return None,None,None,None
    param_names_sigbg, param_values_sigbg, param_errors_sigbg = read_bruout(sigbg_file)
    
    # Now get the bg sideband data
    bg_file="/".join(subdir.split('/')[:-2])+"/outObsBins_sdbnd_bg"
    bg_file+="/"+get_directories_with_values(bg_file, x_values)[0]+"/ResultsHSMinuit2.root"
    if(not os.path.exists(bg_file)):
        return None,None,None,None
    param_names_bg, param_values_bg, param_errors_bg = read_bruout(bg_file)

    # Return
    param_names = param_names_sigbg
    param_values = []
    param_errors = []
    for pn,pv_sigbg,pv_bg,pe_sigbg,pe_bg in zip(param_names,param_values_sigbg,param_values_bg,
                                               param_errors_sigbg,param_errors_bg):
        if("Yld" in pn):
            continue
        
        param_values.append([pv_sigbg[0]/purity-(1-purity)*pv_bg[0]/purity])
        param_errors.append([math.sqrt((pe_sigbg[0]/purity)**2+((1-purity)*pe_bg[0]/purity)**2)])
    
    
    return [round(float(x),2) for x in x_values], param_names, param_values, param_errors

def create_asym_yaml(project_dir="",project_name="",path=""):
    if(path==""):
        path=f"/work/clas12/users/gmat/clas12/clas12_dihadrons/projects/{project_name}/volatile/asym"
    data = {}
    nested_dict = {}
    current_dict = nested_dict
    L=0
    Lmax=100
    for root, dirs, files in os.walk(path):
        for dir in dirs:
            subdir = os.path.join(root,dir)
            split_subdir = subdir.split("/")
            asym_index = split_subdir.index("asym")
            result = split_subdir[asym_index+1:]
            if any(element.startswith("out") for element in result) and \
               result.index(next(filter(lambda x: x.startswith("out"), result))) < len(result) - 1 and \
                ".ipynb_checkpoints" not in result:

                string = ""
                # Now get sPlot or Sdbnd data
                if("outObsBins_splot" in result):
                    x_values, param_names, param_values, param_errors = read_splot(subdir)
                    result[-2]="splot"
                elif("outObsBins_splot_sig" in result):
                    x_values, param_names, param_values, param_errors = read_splot(subdir)
                    result[-2]="splot_sig"
                elif("outObsBins_splot_bg" in result):
                    x_values, param_names, param_values, param_errors = read_splot(subdir)
                    result[-2]="splot_bg"
                elif("outSdbndBins" in result):
                    x_values, param_names, param_values, param_errors = read_sideband(subdir)
                    result[-2]="sideband"
                elif("outObsBins" in result):
                    x_values, param_names, param_values, param_errors = read_splot(subdir)
                    result[-2]="standard"
                else:
                    continue
                data_list=result
                
                if(x_values==None):
                    continue
                
                for i, key in enumerate(data_list[:-1]):
                    if(not key in current_dict):  
                        current_dict[key] = {}
                    if(i==len(data_list[:-1])-1): break
                    current_dict = current_dict[key]

                
                if 'x' in current_dict[key]:
                    current_dict[key]["x"].extend(x_values)
                    for p, v, e in zip(param_names, param_values, param_errors):
                        if p in current_dict[key]:
                            current_dict[key][p]["value"].extend(v)
                            current_dict[key][p]["error"].extend(e)
                        else:
                            current_dict[key][p] = {"value": v, "error": e}
                else:
                    current_dict[key] = {"x": x_values}
                    for p, v, e in zip(param_names, param_values, param_errors):
                        current_dict[key][p] = {"value": v, "error": e}
                current_dict = nested_dict
                
    with open(f"{project_dir}/{project_name}/dihadron_binning.yaml", "w") as file:
        yaml.dump(nested_dict, file)

def get_asym_yaml(project_dir,project_name):
    # Load YAML data from file
    with open(f"{project_dir}/{project_name}/dihadron_binning.yaml", "r") as file:
        data = yaml.load(file, Loader=yaml.FullLoader)
    return data

def get_data_from_yaml(data, headers, par):
    current_dict = data
    for header in headers:
        if(header==""):
            continue
        current_dict = current_dict[header]

    x = current_dict["x"]
    y = current_dict[par]["value"]
    yerr = current_dict[par]["error"]
    
    # Sort the x, y, and yerr arrays by x
    sorted_indices = sorted(range(len(x)), key=lambda i: x[i])
    x = [x[i] for i in sorted_indices]
    y = [y[i] for i in sorted_indices]
    yerr = [yerr[i] for i in sorted_indices]

    return x, y, yerr

def print_dict_keys(dict_obj, level=0):
    # Print keys of dictionary at current level
    key0=""
    for key in dict_obj.keys():
        print("  " * level + key)
        if(key0==""):
            key0=key
    # If value is another dictionary, recursively print its keys
    if isinstance(dict_obj[key0], dict):
        print_dict_keys(dict_obj[key0], level+1)
        
# Get the injection functions for the specific plot
def get_inject_func(binyaml,plot):
    binningStructures = None

    with open(binyaml) as f:
        binningStructures = yaml.safe_load(f)

    sigfunc_dict = {}
    bgfunc_dict = {}

    for structure in binningStructures['binningStructures']:
        name = str(structure['name'])
        sigfuncs = structure.get('inject_sigfuncs',[])
        bgfuncs = structure.get('inject_bgfuncs',[])

        if(name==plot):
            return sigfuncs,bgfuncs
    

def get_inject_plot(func="",xmin=0,xmax=1,y=0): # y is reserved for 2d binning
    x = np.linspace(xmin,xmax,100)
    z = eval(func.replace("sin","np.sin"))
    return np.array(x),np.array(z)


# Returns a dictionary of the LaTeX code for each of the theta integrate modulations
def get_modulations(L):
    char_vec = []
    str_vec = []
    base_str_vec=[]
    for l in range(0, L+1):
        for m in range(1, l+1):
            base_str_vec.append(f"\sin({m}\phi_h-{m}\phi_R)")
            if(m==1):
                string = "\sin(" + "\phi_{h}-" + "\phi_{R})"
            else:
                string = "\sin(" + str(m) +"\phi_{h}-" + str(m) +"\phi_{R})"
            str_vec.append(string)
        for m in range(-l, l+1):
            base_str_vec.append(f"\sin({1-m}\phi_h+{m}\phi_R)")
            if(m==1):
                string = "\sin("+"\phi_{R})"
            elif(m==2):
                string = "\sin("  + "-\phi_{h}+" + str(m) +"\phi_{R})"
            elif(m==0):
                string = "\sin("+"\phi_{h})"
            elif(m==-1):
                string = "\sin(" + str(1-m) + "\phi_{h}-" + "\phi_{R})"
            elif(m<0):
                string = "\sin(" + str(1-m) + "\phi_{h}" + str(m) +"\phi_{R})"
            else:
                string = "\sin(" + str(1-m) + "\phi_{h}+" + str(m) +"\phi_{R})"
            str_vec.append(string)
    # Remove duplicate entries
    idx=np.argsort(base_str_vec)
    
    str_vec=list(np.array(str_vec)[idx])
    str_vec = list(dict.fromkeys(str_vec))
    cidx = 0
    for c in range(ord('A'), ord('A')+len(str_vec)):
        string = ""
        string += chr(c)
        char_vec.append(string)
        cidx += 1
        
    data = (char_vec, str_vec)
    mod_data = {}

    for letter, sin_value in zip(data[0], data[1]):
        mod_data[letter] = sin_value
    return mod_data



# Returns a dictionary of the LaTeX code for each of the theta integrate modulations
# Add the delta phi modulations
def get_2h_modulations(L):
    mod_data = get_modulations(L)
    
    # Find the last key in the dictionary
    last_key = max(mod_data.keys())

    # Get the next uppercase character after the last key
    next_character = chr(ord(last_key) + 1) if last_key != 'Z' else 'A'
    # Get the next next uppercase character after the last key
    nextnext_character = chr(ord(last_key) + 2) if last_key != 'Z' else 'A'
    
    mod_data[next_character] = "\sin(\Delta\phi)"
    mod_data[nextnext_character] = "\sin(2\Delta\phi)"
    
    return mod_data

############
# PLOTTING
############

def get_color_from_channel(channel,colors):
    dihadron_pair=channel[2]
    #colors=["gold","blue","orangered","maroon","green","black"]
    if(dihadron_pair=="piplus_pi0"):
        return colors[0]
    elif(dihadron_pair=="piminus_pi0"):
        return colors[1]
    elif(dihadron_pair=="piplus_piminus"):
        return colors[2]
    elif(dihadron_pair=="piplus_piplus"):
        return colors[3]
    elif(dihadron_pair=="piminus_piminus"):
        return colors[4]
    elif(dihadron_pair=="pi0_pi0"):
        return colors[5]
    else:
        return colors[0]
    
def get_markerstyle_from_channel(channel):
    dihadron_pair=channel[2]
    machine_learning=channel[4]
    if(not "pi0" in dihadron_pair):
        return "o"
    elif(machine_learning=="ML"):
        return "o"
    else:
        return "o"

on_plot_legend = True

def plot_channels(data, channels, drop_edges=False, my_list=None, out_plot=None, extra_legends=None, y_lim=None, custom_legends=None, make_title=True):
    num_plots = 7
    
    # Set up the figure and subplots
    if my_list is None:
        my_list = list(range(num_plots))

    if on_plot_legend:
        figsize = (2*(len(my_list)), 2)
        width_ratios = [(1 if len(my_list)>1 else 1.5)]*len(my_list) 
    else:
        figsize = (2*(len(my_list)+1.5), 2)
        width_ratios = [(1 if len(my_list)>1 else 1.5)]*len(my_list) + [2]
    
    
    fig, axs = plt.subplots(1, len(my_list)+(0 if on_plot_legend else 1), figsize=figsize, dpi=150, sharey=True, sharex=True, gridspec_kw={'wspace': 0.1, 'width_ratios': width_ratios})
    #fig.suptitle('My Plot Title', fontsize=16, y=1.1, x=len(my_list) / (2 * (len(my_list) + 1)))
    
    for j, channel in enumerate(channels):
        channel = channels[j]
        channel_label = format_channel_label(channel,channels)
        extra_legend = ("" if extra_legends==None else extra_legends[j])
        custom_legend = ("" if custom_legends==None else custom_legends[j])
        plot_data(data, axs, channel, channel_label, j, channels, my_list, drop_edges, extra_legend, custom_legend,make_title)
    
    format_subplots(axs, len(my_list))
    if(on_plot_legend):
        create_legend_onplot(axs, len(my_list))
    else:
        create_legend_subplot(axs, len(my_list))
    
    if(y_lim!=None):
        axs[0].set_ylim(y_lim[0],y_lim[1])
        
    if(out_plot!=None):
        plt.savefig(out_plot)
        plt.close()
    
def dup_items(items, channels):
    if(type(items)!=list):
        items=[items]
    
    occurences = 0
    for channel in channels:
        all_found=True
        for item in items:
            if(not item in channel):
                all_found=False
        if(all_found):
            occurences+=1
    print(occurences)
    if(occurences>1):
        return True
    else:
        return False

def format_dihadron_label(dihadron_pair):
    return dihadron_pair.replace("piplus", "$\pi^{+}$").replace("piminus", "$\pi^{-}$").replace("pi0", "$\pi^{0}$").replace("_", "")

def format_channel_label(channel,channels):
    headers = channel
    version = headers[0]
    dihadron_pair = headers[2]
    machine_learning = headers[4]
    
    version = version.replace("Fall2018Spring2019_RGA_inbending","inb. rg-a")
    version = version.replace("Fall2018_RGA_inbending","inb. f18 rg-a").replace("Fall2018_RGA_outbending","outb. f18 rg-a").replace("Spring2019_RGA_inbending","inb. sp19 rg-a").replace("MC_RGA_inbending","inb. MC rg-a").replace("MC_RGA_outbending","outb. MC rg-a")
    version = version.replace("Spring2020_RGB_inbending","inb. sp20 rg-b").replace("Fall2019_RGB_outbending","outb. f19 rg-b").replace("Spring2019_RGB_inbending","inb. sp19 rg-b").replace("MC_RGB_inbending","inb. MC rg-b").replace("MC_RGB_outbending","outb. MC rg-b")
    plot_label = format_dihadron_label(dihadron_pair)
    
    if "pi0" in dihadron_pair:
        plot_label += (" ML" if machine_learning=="ML" else (" Trad." if machine_learning=="noML" else " "))
        
    #if dup_items(dihadron_pair,channels):
    plot_label += f" ({version})"
        
    if "splot_sig" in headers[6]:
        plot_label += " sPlot sig"
    elif "splot_bg" in headers[6]:
        plot_label += " sPlot bkg"
    elif "sideband" in headers[6]:
        plot_label += " sideband sig"
    return plot_label

def plot_data(data, axs, channel, channel_label, channel_index, channels, my_list, drop_edges, extra_legend, custom_legend, make_title):
    dihadron_pair = channel[2]
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color'][:6]
    # Set color for the plot based on the dihadron pair and channel index
    # Color Code Commented on June 29th for the purpose of making the sideband+splot same colors
    #if(dup_items(dihadron_pair,channels)):
    #    color = colors[channel_index]
    #else:
    #    color = get_color_from_channel(channel,colors)
    color = get_color_from_channel(channel,colors)
    # Get the marker style for the plot based on the channel
    markerstyle = get_markerstyle_from_channel(channel)
    
    # Store plot arguments in a dictionary
    plot_kwargs = {
        'fmt': markerstyle,
        'color': color,
        'markersize': 3,
        'capsize': 1.2,
        'label': r"{} {}".format(channel_label, extra_legend)
    }
    
    if channel[6] == "sideband":
        plot_kwargs['markerfacecolor'] = 'white'
    
    if custom_legend:
        plot_kwargs['label'] = format_dihadron_label(custom_legend)
    
    for ix_axis,i in enumerate(my_list):
        # Get data for the plot
        x, y, yerr = get_data_from_yaml(data, channel, get_modulations(2)[0][i])
        # If drop_edges is True, remove the first and last elements of x, y, and yerr
        if drop_edges:
            x, y, yerr = x[1:-1], y[1:-1], yerr[1:-1]

        # Plot the data with the given arguments
        axs[ix_axis].errorbar(x, y, yerr=yerr, **plot_kwargs)
        if make_title:
            axs[ix_axis].set_title(r"${}$".format(get_modulations(2)[1][i]))
        else:
            axs[ix_axis].set_title("")
        
        # If the channel is MC, plot the injection function
        if "MC" in channel[0]:
            sfunc, bfunc = get_inject_func("/work/clas12/users/gmat/clas12/clas12_dihadrons/utils/binning_files/Binning_1d_only.yaml", f"{channel[3]}_binned")
            func = sfunc[i] if ("_sig" in channel[6] or not "pi0" in channel[2] or "sideband" in channel[6]) else bfunc[i]
            x_inject, y_inject = get_inject_plot(func, xmin=np.amin(x), xmax=np.amax(x))
            axs[ix_axis].plot(x_inject, y_inject, color=color, linestyle="dashed", linewidth=1,label=r"{} inject".format(format_dihadron_label(channel[2])))
            
        # Add horizontal dashed line at y=0
        axs[ix_axis].axhline(y=0, color='black', linestyle='--')
        # Turn on the grid
        axs[ix_axis].grid(True)
        # Set the x-axis label based on the channel
        xlabel = channel[3].replace("Mh", "$M_{h}[GeV]$").replace("pTtot","$p_{T}[GeV]$").replace("Mx","$M_{X}[GeV]$")
        axs[ix_axis].set_xlabel(r"{}".format(xlabel))
        
def format_subplots(axs, num_plots):
    axs[0].set_ylabel(r"$A_{LU}$")
    ymax = np.amax(np.abs(axs[0].get_ylim()))
    axs[0].set_ylim(-ymax, ymax)
    axs[0].set_xlim(axs[0].get_xlim()[0]-0.1,axs[0].get_xlim()[1]+0.1)
    for i in range(1, num_plots+(0 if on_plot_legend else 1)):
        axs[i].tick_params(axis='y', which='both', labelleft=False)

def create_legend_subplot(axs, num_plots):
    handles, labels = axs[0].get_legend_handles_labels()

    legend = axs[-1].legend(handles, labels, loc='upper center', ncol=1, bbox_to_anchor=(0.5, -0.15),fontsize=12)
    legend.set_bbox_to_anchor((0.6, 0.9))
    axs[-1].text(0.5, 0.92, 'Key', ha='center', transform=axs[-1].transAxes)
    axs[-1].axis('off')
    axs[-1].grid(False)

    for i in range(0, num_plots):
        axs[i].legend().remove()
        
def create_legend_onplot(axs, num_plots):
    
    handles, labels = axs[0].get_legend_handles_labels()
    labels = [label.split(" ")[0] for label in labels]
    
    # Remove duplicates
    seen_labels = set()
    unique_labels = []
    unique_handles = []

    for label, handle in zip(labels, handles):
        if label not in seen_labels:
            unique_labels.append(label)
            unique_handles.append(handle)
            seen_labels.add(label)
    
    legend = axs[0].legend(unique_handles, unique_labels, loc='upper left', ncol=1, fontsize=10,frameon=True,bbox_to_anchor=(0.03, 0.97),
                          borderpad=0.15)
    legend.get_frame().set_facecolor('white')  # Set the legend background color
    legend.get_frame().set_edgecolor('grey')  # Set the legend border color
    legend.get_frame().set_linewidth(1.5)  # Set the border thickness
    legend.get_frame().set_boxstyle('round,pad=0.5')  # Set the box style with rounded corners
