import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import argparse
import os

import amino_acid_analysis.analysis as analysis
import amino_acid_analysis.find_correlation as correlation
from amino_acid_analysis.properties import hydrophobic_residues, hbond_residues

def plot_all(args):
    """Make all the plots"""
    
    cwd = os.getcwd() #Assumes you are in the temperature directory

    params_list = [] #List of all params

    for i in args.iter:
        #Find and load files
        os.chdir("%s/iteration_%s/"%(cwd, args.iter[i], args.temp)) ##Hardcoded directory search line, will generalize at some point
        atom_index, res_index, atom_type, res_id = analysis.read_pdb("%s"%args.atom_file)
        pairs = np.loadtxt("%s"%args.pairs_file, usecols=(0, 1))
        params = np.loadtxt("%s"%args.params_file)
        params_list.append(params)
        
        os.chdir(cwd)

        if not os.path.isdir("contact_plots"):
            os.mkdir("contact_plots")
    
        os.chdir("contact_plots")
        
        CACB_magnitude(pairs, params, res_index, atom_type, iter=i)
        polarity_magnitude(pairs, params, res_index, atom_type, res_id, iter=i)
        CB_contact_type(pairs, params, res_index, atom_type, res_id, iter=i)
        
    #Find correlated groups, plot against CA/CB and polarity
    pairs_corr, params_corr = correlation.generate_correlated_groups(pairs, params_list)
    ##Set as iter+1 temporarily, generalize plotting functions at some point
    CACB_magnitude(pairs_corr, params_corr, res_index, atom_type, iter=max(args.iter)+1)
    polarity_magnitude(pairs_corr, params_corr, res_index, atom_type, res_id, iter=max(args.iter)+1)
    
###Convenience functions to generate specific types of plots###
def CACB_magnitude(pairs, params, res_index, atom_type, iter=0):
    """Plot magnitude of CA epsilons in upper triangle, magnitude of CB epsilons in lower triangle"""
    
    ca_pairs, ca_params, cb_pairs, cb_params = analysis.separate_cacb(pairs, params, atom_type)
    ca_res_index = analysis.get_residue_indices(ca_pairs, res_index)
    cb_res_index = analysis.get_residue_indices(cb_pairs, res_index)
    
    plot_color(ca_res_index, cb_res_index, ca_params, cb_params, title="Iteration %.0f"%iter, savename="epsilons_iter%.0f"%iter)
    plot_spread((ca_params, cb_params), ("CA", "CB"), title="Epsilon Spread", savename="cacb_eps_spread_iter_%.0f.png"%iter)

def polarity_magnitude(pairs, params, res_index, atom_type, res_id, iter=0):
    """Plot magnitude of hydrophobic epsilons in upper triangle, magnitude of hydrogen bonding epsilons in lower triangle"""
    
    ca_pairs, ca_params, cb_pairs, cb_params = analysis.separate_cacb(pairs, params, atom_type)

    hydrophobic_pairs, hydrophobic_params, hbond_pairs, hbond_params = analysis.separate_polarity(cb_pairs, cb_params, res_id=res_id)
    hydrophobic_res_index = analysis.get_residue_indices(hydrophobic_pairs, res_index)
    hbond_res_index = analysis.get_residue_indices(hbond_pairs, res_index)

    plot_color(hydrophobic_res_index, hbond_res_index, hydrophobic_params, hbond_params, title="Iteration %.0f"%iter, savename="polarity_epsilons_iter%.0f"%iter)
    plot_spread((hydrophobic_params, hbond_params), ("Hydrophobic", "Hydrogen Bond"), title="Epsilon Spread", savename="polarity_eps_spread_iter_%.0f.png"%iter)

def native_magnitude(pairs, params, native_list, iter=0):
    """Plot native contact epsilons in upper triangle, nonnative contact epsilons in lower triangle"""
    native_pairs, native_params, nonnative_pairs, nonnative_params = analysis.separate_native(pairs, params, native_list)
    plot_color(native_pairs, nonnative_pairs, native_params, nonnative_params, title="Iteration %.0f"%iter, savename="native_epsilons_iter_%.0f.png"%iter)
    plot_spread((native_params, nonnative_params), ("native","nonnative"), title="Epsilon Spread", savename="native_eps_spread_iter_%.0f.png"%iter)

def CB_contact_type(pairs, params, res_index, atom_type, res_id, iter=0):
    """Plot magnitude of CB epsilons in upper triangle, contact type in lower triangle"""
    
    ca_pairs, ca_params, cb_pairs, cb_params = analysis.separate_cacb(pairs, params, atom_type)
    cb_res_index = analysis.get_residue_indices(cb_pairs, res_index)
    contact_type = analysis.sort_polarity(cb_pairs, res_id)
    
    plot_color(cb_res_index, cb_res_index, cb_params, contact_type, title="Iteration %.0f"%iter, savename="contact_types_iter%.0f"%iter)

###General plotting functions
def plot_color(pairs_U, pairs_L, params_U, params_L, title="Contacts", savename="contacts.png"):
    """General function to generate contact map plot"""

    x_U = pairs_U[:,0]
    y_U = pairs_U[:,1]
    z_U = params_U

    x_L = pairs_L[:,0]
    y_L = pairs_L[:,1]
    z_L = params_L

    # Set constant zmin, zmax
    zmin = 0.1
    zmax = 2.0
    params_U[params_U < zmin] = zmin
    params_L[params_L < zmin] = zmin
    params_U[params_U > zmax] = zmax
    params_L[params_L > zmax] = zmax

    plt.figure()
    
    cp = plt.scatter(y_L, x_L, s=6, c=z_L, marker='o', linewidth=0., vmin=zmin, vmax=zmax)
    cp = plt.scatter(x_U, y_U, s=6, c=z_U, marker='o', linewidth=0., vmin=zmin, vmax=zmax)
    cb = plt.colorbar(cp)

    max_U = np.max(x_U)
    max_L = np.max(x_L)
    if max_U > max_L:
        maxval = max_U
    else:
        maxval = max_L
    plt.axis([0, maxval, 0, maxval])
    plt.xlabel("Residue i")
    plt.ylabel("Residue j")
    plt.title("%s"%title)

    plt.savefig("%s"%savename)
    plt.close()

def plot_spread(params_list, labels_list, title="Epsilon Spread", savename="eps_spread.png"):
    """General function to plot and compare spreads"""
    plt.figure()
    for i in range(len(params_list)):
        plt.hist(params_list[i], bins=50, alpha=0.5, label="%s"%labels_list[i])
    plt.xlabel("Epsilon Value")
    plt.ylabel("Count")
    plt.title("%s"%title)
    plt.legend()
    plt.savefig("%s"%savename)

def get_args():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--atom_file", default="Native.pdb", type=str, help="File containing atom data")
    parser.add_argument("--pairs_file", default="pairwise_params", type=str, help="File containing contact pairs by atom index")
    parser.add_argument("--params_file", default="model_params", type=str, help="File containing contact parameters")
    parser.add_argument("--temp", type=int, help="temperature directory")
    parser.add_argument("--iter", nargs="+", type=int, help="iteration")
    
    args = parser.parse_args()
    
    return args

if __name__=="__main__":
   
    args = get_args()
    plot_all(args)
