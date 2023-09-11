#!/usr/bin/env python
# coding: utf-8

"""
Script to estimate interface pLDDT and guide on use of docking
protocols in the following step of the pipeline.

Author: Ameya Harmalkar (aharmal1@jhu.edu)
"""

import pandas as pd
import numpy as np
import argparse
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.Superimposer import Superimposer
import urllib.request

def get_interface_residues(pdb_file, partners="A_B", cutoff=10.0):
    """Obtain the interface residues in a particular predicted structure. 
    Args
        pdb_file (str): Path to a PDB file
        partners (str): String annotating the docking partners. Partners are 
                        defined as <partner chain 1>_<partner chain 2>. for e.g.
                        A_B, A_HL, ABC_DE, and so on. Note this is applicable for
                        proteins and peptides with canonical amino acids only.
        cutoff (float): Cut-off to determine the interface residues. Default is 
                        10 Angstorms.
    """
    parser = PDBParser()
    structure = parser.get_structure("model", pdb_file)

    partner1 = partners.split("_")[0]
    partner2 = partners.split("_")[1]

    interface_residues = {'partner1': [], 'partner2': []}

    for chain1 in partner1:
        for chain2 in partner2:
            chain1_residues = set(structure[0][chain1].get_residues())
            chain2_residues = set(structure[0][chain2].get_residues())
            int1_residues = []
            int2_residues = []

            for residue1 in chain1_residues:
                for residue2 in chain2_residues:
                    distance = min(
                        [np.linalg.norm(atom1.get_coord() - atom2.get_coord())
                         for atom1 in residue1.get_atoms()
                         for atom2 in residue2.get_atoms()])

                    if distance <= cutoff:
                        int1_residues.append((chain1,residue1))
                        int2_residues.append((chain2,residue2))

            interface_residues['partner1'].extend(list(set(int1_residues)))
            interface_residues['partner2'].extend(list(set(int2_residues)))

    interface_residues['partner1'] = list(set(interface_residues['partner1']))
    interface_residues['partner2'] = list(set(interface_residues['partner2']))

    return interface_residues


def get_interface_residue_b_factors(pdb_file, partners="A_B", cutoff = 10):
    """
    Returns a dictionary of interface residues and their B-factor values.
    :param pdb_file: The path to the PDB file.
    :param chain1: The ID of the first chain.
    :param chain2: The ID of the second chain.
    :param cutoff_distance: The cutoff distance for identifying interface residues.
    """
    structure = PDBParser().get_structure('pdb', pdb_file)
    model = structure[0]
    b_factors = {}
    
    interface_residues = get_interface_residues(pdb_file, partners, cutoff )
    
    for k, v in interface_residues.items():
        for residue in v:
            for atom in residue[-1]:
                if ( atom.get_name() == 'CA' ):
                    b_factors[(k, residue[0], residue[-1].get_id()[1], atom.get_name())] = atom.get_bfactor()
            
    return b_factors


def get_b_factors(pdb_file):
    """
    Returns the average pLDDT (or b-factor) summed up across all the residues.
    Note that in AlphaFold predictions, the b-factor column holds the pLDDT
    information.
    Args
        pdb_file (str) : The path to the PDB file.
    Returns
        avg_b_factor (float) : Average value of B-factor/pLDDT.
    """
    structure = PDBParser().get_structure('pdb', pdb_file)
    model = structure[0]
    model_residues = model.get_residues()
    b_factors = {}
    
    for residue in model_residues:
        for atom in residue:
            if ( atom.get_name() == 'CA' ):
                b_factors[( residue.get_id()[1], atom.get_name())] = atom.get_bfactor()
                
    avg_b_factor = np.mean(list(b_factors.values()))
            
    return avg_b_factor


def structural_confidence(avg_b_factor):
    """Output structural confidence of AlphaFold predictions.
    Args
        avg_b_factor (float) : Average pLDDT value
    """
    if avg_b_factor >= 80:
        print( "STATUS: AlphaFold has HIGH overall confidence in the structure" )
    elif 60 < avg_b_factor < 80:
        print( "STATUS: AlphaFold has LOW overall confidence in the structure. Please check for disordered regions before initiating docking." )
    else:
        print( "STATUS: AlphaFold prediction is possibly incorrect. Docking results would be incorrect." )
    return None


def docking_confidence(interface_plddt):
    """Output docking confidence status for users. 
    Args
        interface_plddt (float) : Interface pLDDT value
    """
    if interface_plddt >= 85:
        print("STATUS: Docking site is predicted with higher confidence. Use LOCAL Docking and Refinement.")
    else:
        print("STATUS: Docking site has lower confidence. Use GLOBAL Docking and refine using LOCAL docking.")
    return None


def main():
    
    parser = argparse.ArgumentParser(description="Estimate interface pLLDT from AF2 generated structures")
    parser.add_argument( 'model', metavar='<model', type=str, nargs=1, help='Path to model file')
    parser.add_argument( '-partners', type=str, default='A_B', help='option specifying the chainbreak, i.e. protein partners' )
    parser.add_argument( '-cutoff', type=float, default=10, help='Distance cut-off to determine the interface' )
    
    args = parser.parse_args()
    model = args.model[0]
    partners = args.partners
    cutoff = args.cutoff
    
    b_factors = get_interface_residue_b_factors( model, partners, cutoff )
    interface_bfactor = np.mean(list(b_factors.values()))
    net_b_factor = get_b_factors( model )

    print(("Avg pLDDT {:.3f}".format(net_b_factor)))
    structural_confidence(net_b_factor)
    print(("Interface pLDDT {:.3f}".format(interface_bfactor)))
    docking_confidence(interface_bfactor)
    
    
if __name__ == '__main__':
    main() 