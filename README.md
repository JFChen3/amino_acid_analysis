# Amino Acid Analysis

## Purpose: 
Plot and compare epsilons of different types of residue contacts. 

## Description: 
####analysis.py
Functions to get data (atom indices, residue indices, atom types, residue types) from a pdb file and sort contacts. Currently able to identify CA-CA, CB-CB, hydrophobic, and hydrogen bonding interactions.

####plot_contacts.py
General plotting functions for contact map and distribution plots as well as specific functions to generate:
  
  - CA vs CB contact map and epsilon spread
    
  - Hydrophobic vs hydrogen bonding contact map and epsilon spread
    
  - CB contact map vs hydrophobicity

####properties.py
Lists of hydrophobic and hydrogen bonding residues.
