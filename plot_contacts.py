import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import argparse
import os

import amino_acid_analysis.analysis as analysis
from amino_acid_analysis.properties import hydrophobic_residues, hbond_residues

def plot_all(args):
    """Make all the plots"""
    
    cwd = os.getcwd() #Assumes you are in the temperature directory

    for i in args.iter:
        #Find and load files
        os.chdir("%s/1PB7/iteration_%s/%s_0"%(cwd, args.iter[i], args.temp))
        atom_index, res_index, atom_type, res_id = analysis.read_pdb("%s"%args.atom_file)
        pairs = np.loadtxt("%s"%args.pairs_file, usecols=(0, 1))
        params = np.loadtxt("%s"%args.params_file)

        os.chdir(cwd)

        if not os.path.isdir("contact_plots"):
            os.mkdir("contact_plots")
    
        os.chdir("contact_plots")
        
        CACB_magnitude(pairs, params, res_index, atom_type, iter=i)
        polarity_magnitude(pairs, params, res_index, res_id, iter=i)
        CB_contact_type(pairs, params, res_index, atom_type, res_id, iter=i)

###Convenience functions to generate specific types of plots
def CACB_magnitude(pairs, params, res_index, atom_type, iter=0):
    """Plot magnitude of CA epsilons in upper triangle, magnitude of CB epsilons in lower triangle"""
    
    ca_pairs, ca_params, cb_pairs, cb_params = analysis.separate_cacb(pairs, params, atom_type)
    ca_res_index = analysis.get_residue_indices(ca_pairs, res_index)
    cb_res_index = analysis.get_residue_indices(cb_pairs, res_index)
    
    analysis.plot_color(ca_res_index, cb_res_index, ca_params, cb_params, title="Iteration %.0f"%iter, savename="epsilons_iter%.0f"%iter)
    
    #Plot spreads
    analysis.plot_spread((ca_params, cb_params), ("CA", "CB"), title="Epsilon Spread", savename="cacb_eps_spread_iter_%.0f.png"%iter)

def polarity_magnitude(pairs, params, res_index, res_id, iter=0):
    """Plot magnitude of hydrophobic epsilons in uppper triangle, magnitude of hydrogen bonding epsilons in lower triangle"""
    
    hydrophobic_pairs, hydrophobic_params, hbond_pairs, hbond_params = analysis.separate_polarity(pairs, params, res_id=res_id)
    hydrophobic_res_index = analysis.get_residue_indices(hydrophobic_pairs, res_index)
    hbond_res_index = analysis.get_residue_indices(hbond_pairs, res_index)

    analysis.plot_color(hydrophobic_res_index, hbond_res_index, hydrophobic_params, hbond_params, title="Iteration %.0f"%iter, savename="epsilons_iter%.0f"%iter)

    #Plot spreads
    analysis.plot_spread((hydrophobic_params, hbond_params), ("Hydrophobic", "Hydrogen Bond"), title="Epsilon Spread", savename="polarity_eps_spread_iter_%.0f.png"%iter)

def CB_contact_type(pairs, params, res_index, atom_type, res_id, iter=0):
    """Plot magnitude of CB epsilons in upper triangle, contact type in lower triangle"""
    
    ca_pairs, ca_params, cb_pairs, cb_params = analysis.separate_cacb(pairs, params, atom_type)
    cb_res_index = analysis.get_residue_indices(cb_pairs, res_index)
    contact_type = analysis.sort_polarity(cb_pairs, res_id)
    
    analysis.plot_color(cb_res_index, cb_res_index, cb_params, contact_type, title="Iteration %.0f"%iter, savename="contact_types_iter%.0f"%iter)

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