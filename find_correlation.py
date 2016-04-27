import numpy as np

def generate_correlated_groups(contact_pairs, params_list):
    #Walks through steps
    
    params_matrix = construct_params_matrix(params_list)
    cov_matrix = calc_covariance_matrix(params_matrix)
    groups = group_contacts(cov_matrix, 0.9) ###Arbitrarily set Rmin to 0.9, set as argument later
    pairs_plot, params_plot = get_plottable_contacts(contact_pairs, groups)

    return pairs_plot, params_plot

def construct_params_matrix(params_list):
    #Construct parameter matrix
    matrix = np.zeros((np.size(params_list[0]), np.size(params_list[0])))
    for i in range(len(params_list)):
        matrix[:,i] = params_list[i]

    return matrix

def calc_covariance_matrix(params_matrix):
    #Calculate normalized covariance matrix
        
    cov_matrix_normalized = np.corrcoef(matrix)
    
    return cov_matrix_normalized

def group_contacts(cov_matrix, R_min):
    #Group together contacts that change similarly, returns indices of correlated contacts

    n_pairs = np.shape(cov_matrix)[0]

    correlated_pairs = np.where(cov_matrix > R_min)
    
    current_group = 0 #Track number of groups
    
    groups = [[0]]

    for i in range(np.size(correlated_pairs[0])): #Iterate through correlated_pairs
        if correlated_pairs[0][i] < correlated_pairs[1][i]: #Exclude pairs where second index is larger than first to not double count
            if correlated_pairs[0][i] == correlated_pairs [0][i-1]: #If first entry is still the same, add to same group
                groups[current_group].append(correlated_pairs[1][i])            
            else: #if different first entry, make new group
                groups.append([])
                current_group +=1
                groups[current_group].append(correlated_pairs[0][i])
                groups[current_group].append(correlated_pairs[1][i])

    return groups
    
def get_plottable_contacts(pairs, groups):
    #pairs are actual residue indices
    #returns groups with color coded params
    
    pairs_use = np.zeros((0,2))
    pair_markers = []
    n_groups = len(groups)
    spacing = 2.0/(n_groups - 1.0)
    
    for i in range(n_groups):
        for j in range(len(groups[i])):
            pairs_use=np.vstack(pairs_use, pairs[groups[i][j]])
            pair_markers.append(i*spacing)
    
    pair_markers = np.array(pair_markers)
    
    return pairs_use, pair_markers