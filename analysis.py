import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from properties import hydrophobic_residues, hbond_residues

def read_pdb(filename):
    """Get atom indices, residue indices, atom types, and residue types from Native.pdb file"""

    info_pdb = open("%s"%filename, "r")

    atom_index = []
    res_index = []
    atom_type = []
    res_id = []

    for line in info_pdb: #Iterate through lines of pdb file
        info = line.strip().split()
        if info[0] == "END":
            break
        atom_index.append(info[1])
        res_index.append(info[5])
        atom_type.append(info[2])
        res_id.append(info[3])

    return atom_index, res_index, atom_type, res_id

def get_residue_indices(atom_pairs, res_index):
    """Gets residue indices from atom indices"""
    
    residue_pairs = np.zeros((np.shape(atom_pairs)))
    
    for i in range(np.shape(atom_pairs)[0]):
        for j in range(np.shape(atom_pairs)[1]):
            residue_pairs[i][j] = res_index[int(atom_pairs[i][j])-1]

    return residue_pairs

def sort_polarity(pairs, res_id):
    """Sort contacts into hydrophobic, dipole, and hydrogen bond interactions.
       Assigns a 0 for hydrophobic, 2 for hydrogen bonds, and 1 for everything else
       Inputs - pairs: array of atom pairs
                res_id: 3-letter residue IDs listed by atom index"""

    atom_i = pairs[:,0] #Assumes pairs file gives atom indices
    atom_j = pairs[:,1]

    n_pairs = np.shape(pairs)[0]
    
    contact_type = np.ones(n_pairs)

    for i in range(n_pairs):
        if res_id[int(atom_i[i])-1] in hydrophobic_residues and res_id[int(atom_j[i])-1] in hydrophobic_residues:
            contact_type[i] -= 1 #Set hydrophobic interactions to 0
        if res_id[int(atom_i[i])-1] in hbond_residues and res_id[int(atom_j[i])-1] in hbond_residues:
            contact_type[i] += 1 #Set hydrogen bond interactions to 2
        else:
            pass #Set everything else to 1

    contact_type = np.array(contact_type)
        
    return contact_type

def separate_cacb(pairs, params, atom_type):
    """Separate contacts into CACA and CBCB"""

    ca_pairs = []
    ca_params = []
    cb_pairs = []
    cb_params = []

    #Count CACA and CBCB interactions
    ca_count = 0
    cb_count = 0

    for i in range(np.shape(pairs)[0]):
        if not atom_type[int(pairs[i,0])-1] == atom_type[int(pairs[i,1])-1]:
            print "Warning: Found CA-CB contact between atoms %.0f and %.0f"%(pairs[i,0], pairs[i,1])

        if atom_type[int(pairs[i,0])] == 'CA':
            ca_pairs.append(pairs[i,:])
            ca_params.append(params[i])
            ca_count += 1
        elif atom_type[int(pairs[i,0])] == 'CB':
            cb_pairs.append(pairs[i,:])
            cb_params.append(params[i])
            cb_count += 1
        else:
            raise ValueError("Atom not CA or CB")

    ca_pairs = np.reshape(ca_pairs, (ca_count,2)) #Reshape into nx2 matrices
    cb_pairs = np.reshape(cb_pairs, (cb_count,2))

    return ca_pairs, ca_params, cb_pairs, cb_params

def separate_polarity(pairs, params, contact_type=None, res_id=None):
    """Separate contacts into hydrophobic and hydrogen bond interactions"""
        
    if contact_type == None:
        if res_id == None:
            raise IOError("Must provide contact_type or res_id")
        else:
            contact_type = sort_polarity(pairs, res_id)
   
    hydrophobic_pairs = []
    hydrophobic_params = []
    hbond_pairs = []
    hbond_params = []

    #Count hydrophobic and hydrogen bond interactions
    hydrophobic_count = 0
    hbond_count = 0

    for i in range(np.shape(pairs)[0]):

        if contact_type[i] == 0:
            hydrophobic_pairs.append(pairs[i,:])
            hydrophobic_params.append(params[i])
            hydrophobic_count += 1
        elif contact_type[i] == 2:
            hbond_pairs.append(pairs[i,:])
            hbond_params.append(params[i])
            hbond_count += 1
        else:
            print "Contact between %.0f and %.0f not characterized as hydrophobic or hydrogen bond"%(pairs[i,0], pairs[i,1])

    hydrophobic_pairs = np.reshape(hydrophobic_pairs, (hydrophobic_count,2)) #Reshape into nx2 matrices
    hbond_pairs = np.reshape(hbond_pairs, (hbond_count,2))
    
    return hydrophobic_pairs, hydrophobic_params, hbond_pairs, hbond_params        
    
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

    maxval = np.max(x_U)
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