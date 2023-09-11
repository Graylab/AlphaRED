
"""
Description: Script to calculate the quantitative metrics such as interface pLDDT,
             avg-pLDDT, Interface contacts, Interface residues.

Author: Ameya Harmalkar (aharmal1@jhu.edu)
"""

import pandas as pd
import numpy as np
import argparse
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.Superimposer import Superimposer
import sys


def get_residues(pdb_file, partners="B_C"):
    """Obtains the residues of the two partners in a PDB file specified by the partners string.
    **Note**: This function assumes that the PDB file contains only one model.
    Args
        pdb_file: str
            The path to the PDB file.
        partners: str
            A string of the form "A_B" where A and B are the chain IDs of the two partners.
    Returns
        partner1_residues: set
            A set of Bio.PDB.Residue.Residue objects for partner 1.
        partner2_residues: set
            A set of Bio.PDB.Residue.Residue objects for partner 2.
    """
    parser = PDBParser()
    structure = parser.get_structure("model", pdb_file)

    partner1 = partners.split('_')[0]
    partner2 = partners.split('_')[1]

    partner1_residues = set()
    partner2_residues = set()
    for i in range(len(partner1)):
        partner1_residues.update( set( structure[0][partner1[i]].get_residues() ) )
    for i in range(len(partner2)):
        partner2_residues.update( set( structure[0][partner2[i]].get_residues() ) )

    if len(partner1_residues) == 0:
        sys.exit("Error: Partner 1 residue set is empty")
    if len(partner2_residues) == 0:
        sys.exit("Error: Partner 2 residue set is empty")

    return partner1_residues, partner2_residues


def get_interface_residues(pdb_file, partners="A_B", cutoff=8):

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


def get_interface_contacts( partner1_residues, partner2_residues, tol = 8.0 ):
    """Returns the interfacial contacts between two partners specified by the residue sets.
    Args
        partner1_residues : set
            Set of Bio.PDB.Residue.Residue objects for partner 1.
        partner2_residues : set
            Set of Bio.PDB.Residue.Residue objects for partner 2.
        tol : float
            Distance cutoff for defining an interface contact.
    Returns
        if_contacts : list
            List of tuples of the form (residue1, residue2) where residue1 and residue2 are
            res IDs obtained from Bio.PDB.Residue.Residue objects.
    """
    if_contacts = []
    if_partner1 = []
    if_partner2 = []

    for residue1 in partner1_residues:
        for residue2 in partner2_residues:
            if residue1.get_resname() == 'GLY' and residue2.get_resname() == 'GLY':
                distance = [np.linalg.norm(atom1.get_coord() - atom2.get_coord())
                        for atom1 in residue1.get_atoms() if atom1.get_name() == 'CA'
                        for atom2 in residue2.get_atoms() if atom2.get_name() == 'CA'][0]
            elif residue1.get_resname() == 'GLY':
                distance = [np.linalg.norm(atom1.get_coord() - atom2.get_coord()) 
                        for atom1 in residue1.get_atoms() if atom1.get_name() == 'CA'
                        for atom2 in residue2.get_atoms() if atom2.get_name() == 'CB'][0]
            elif residue2.get_resname() == 'GLY':
                distance = [np.linalg.norm(atom1.get_coord() - atom2.get_coord()) 
                        for atom1 in residue1.get_atoms() if atom1.get_name() == 'CB'
                        for atom2 in residue2.get_atoms() if atom2.get_name() == 'CA'][0]
            else:
                distance = [np.linalg.norm(atom1.get_coord() - atom2.get_coord()) 
                        for atom1 in residue1.get_atoms() if atom1.get_name() == 'CB'
                        for atom2 in residue2.get_atoms() if atom2.get_name() == 'CB'][0]

            if distance <= tol:
                if_contacts.append((residue1.get_id()[1], residue2.get_id()[1]))
                if_partner1.append((str(residue1.get_id()[1]) + '_' + residue1.get_resname()))
                if_partner2.append((str(residue2.get_id()[1]) + '_' + residue2.get_resname()))
    
    if_part1 = np.unique( if_partner1 )
    if_part2 = np.unique( if_partner2 )

    return len(if_contacts), len(if_part1), len(if_part2)



def get_bfactor(residue):
    """Returns the B-factor/pLDDT of a residue.
    Args
        residue: Bio.PDB.Residue.Residue
            The residue for which the B-factor/pLDDT is to be obtained.
    Returns
        bfactor: float
            The B-factor/pLDDT of the residue.
    """
    for atom in residue.get_atoms():
        if ( atom.get_name() == 'CA' ):
            bfactor = atom.get_bfactor()
    if ( bfactor == None ):
        sys.exit("Error: No CA atom found for residue {}".format(residue.get_id()))
    return bfactor


def get_avg_plddt(pdb_file):
    """
    Returns a dictionary of interface residues and their B-factor values.
    :param pdb_file: The path to the PDB file.
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


def get_interface_residue_plddt( partner1_residues, partner2_residues, tol = 8.0 ):
    """Obtain the pLDDT of the interface between two partners.
    Args
        partner1_residues: set
            Set of Bio.PDB.Residue.Residue objects for partner 1.
        partner2_residues: set
            Set of Bio.PDB.Residue.Residue objects for partner 2.
        tol: float
            Distance cutoff for defining an interface contact (default: 8.0)
    Returns
        partner1_lddt: float
            The pLDDT/bfactor of partner 1.
        partner2_lddt: float
            The pLDDT/bfactor of partner 2.
        interface_lddt: float
            The pLDDT/bfactor of the interface.
    """
    partner1_bfactor = {}
    partner2_bfactor = {}
    for residue1 in partner1_residues:
        for residue2 in partner2_residues:

            distance = min([np.linalg.norm(atom1.get_coord() - atom2.get_coord())
                            for atom1 in residue1.get_atoms()
                            for atom2 in residue2.get_atoms()])
            
            if distance <= tol:
                partner1_bfactor[(residue1.get_id()[1],residue1.get_resname())] = get_bfactor(residue1)
                partner2_bfactor[(residue2.get_id()[1],residue2.get_resname())] = get_bfactor(residue2)

    partner1_lddt = np.mean(list(partner1_bfactor.values()))
    partner2_lddt = np.mean(list(partner2_bfactor.values()))
    interface_lddt = np.mean( list(partner1_bfactor.values()) + list(partner2_bfactor.values())) 
    return partner1_lddt, partner2_lddt, interface_lddt


def get_interface_plddt(pdb_file, partners="A_B", cutoff = 10):
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
    
    interface_bfactor = np.mean(list(b_factors.values()))

    return interface_bfactor


def main():
    parser = argparse.ArgumentParser(description="Estimate interface pLLDT from AF2 generated structures")
    parser.add_argument( 'model', metavar='<model>', type=str, nargs=1, help='Path to model file')
    parser.add_argument( '-partners', type=str, help='Path to partners file (oftenly a Rosetta input)' )
    parser.add_argument( '-part_str', type=str, help='Partners string, e.g. BC_D' )
    parser.add_argument( '-cutoff', type=float, default=8, help='Distance cut-off to determine the interface' )
    
    args = parser.parse_args()
    model = args.model[0]
    cutoff = args.cutoff
    
    # read the text file and get the partners
    partner = ''
    if args.part_str:
        partner = args.part_str
    elif args.partners:
        partner = open(args.partners,'r').readlines()[0].split()[-1]
    
    chain1_residues, chain2_residues = get_residues(model, partner)
    net_res = len(chain1_residues) + len(chain2_residues)
    #partner1_lddt, partner2_lddt, interface_lddt = get_interface_residue_plddt( chain1_residues, chain2_residues, cutoff )
    interface_plddt = get_interface_plddt( model, partner, cutoff )
    if_contacts, if_part1, if_part2 = get_interface_contacts( chain1_residues, chain2_residues, cutoff )
    avg_plddt = get_avg_plddt( model )

    print(("Interface-pLDDT\tAvg-pLDDT\tInterface-contacts\tIF_p1\tIF_p2\tNet_res"))
    print(("{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}".format(interface_plddt, avg_plddt, if_contacts, if_part1, if_part2, net_res)))


if __name__ == '__main__':
    main() 