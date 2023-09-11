#!/usr/bin/env python
# coding: utf-8

"""
Description: Script to calculate RMSD and LDDT values of a model with respect to
             a reference structure

Author: Ameya Harmalkar (aharmal1@jhu.edu)
"""

import pandas as pd
import numpy as np
import pyrosetta
pyrosetta.init('-ignore_unrecognized_res -ignore_zero_occupancy False -mute all')
from pyrosetta import *
from pyrosetta.rosetta.core.scoring import *

from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.Superimposer import Superimposer



def residue_rmsd(model, reference, output=True, outfile='default'):
    pose = pose_from_pdb(model)
    native = pose_from_pdb(reference)
    print(pose.size(),native.size())
    # Define PerResidueRMSDMetric
    perResRMS = pyrosetta.rosetta.core.simple_metrics.per_residue_metrics.PerResidueRMSDMetric()
    perResRMS.set_comparison_pose(native)

    # Define residue selectors
    # reference
    reference = pyrosetta.rosetta.core.select.residue_selector.TrueResidueSelector()
    reference.apply(native)
    # model
    model = pyrosetta.rosetta.core.select.residue_selector.TrueResidueSelector()
    model.apply(pose)

    perResRMS.set_residue_selector(model)
    perResRMS.set_residue_selector_reference(reference)
    perResRMS.set_residue_selector_super(reference)
    perResRMS.set_rmsd_type(rmsd_atoms.rmsd_protein_bb_heavy)
    perResRMS.set_run_superimpose(True)

    df = pd.DataFrame.from_dict(dict(perResRMS.calculate(pose)), orient='index', columns=['rmsd'])
    df['residue'] = df.index
    if output:
        df.to_csv(outfile+'.csv', index=False)
    return df


def residue_lddt(model, reference, output=False, outfile='default'):
    
    pose = pose_from_pdb(model)
    native = pose_from_pdb(reference)
    lddt_calc = pyrosetta.rosetta.core.scoring.lDDT_Calculator()
    lddt_calc.R0(10.0)
    lddt = list(lddt_calc.residue_lDDT(pose, native))
    pdb_res = list(range(1, pose.size() + 1))
    df = pd.DataFrame({'residue': pdb_res, 'lDDT': lddt})
    if output:
        df.to_csv(outfile+'.csv', index=False)
    return df


def af2_plddt(af2_model, partners, output=False, outfile='af2_plddt'):

    parser = PDBParser()
    model = parser.get_structure("model", af2_model)[0]

    partner1 = partners.split("_")[0]
    partner2 = partners.split("_")[1]

    partner1_bfactor = []
    partner2_bfactor = []

    for chains in model.get_chains():
        partner1_bf = []
        partner2_bf = []
        for a in chains:
            if chains.id in partner1:
                ## Chain ID, Residue ID, B-factor
                partner1_bf = [[a.get_full_id()[2], a.get_full_id()[3][1], a.get_bfactor()] for a in chains.get_atoms() if a.name == "CA"]
            elif chains.id in partner2:
                partner2_bf = [[a.get_full_id()[2], a.get_full_id()[3][1], a.get_bfactor()] for a in chains.get_atoms() if a.name == "CA"]
        if partner1_bf != []:
            partner1_bfactor.append(partner1_bf)
        if partner2_bf != []:
            partner2_bfactor.append(partner2_bf)

    p1_dict = { n:[n+1, v[-1]] for n,v in enumerate([item for sublist in partner1_bfactor for item in sublist]) }
    p2_dict = { n:[n+1, v[-1]] for n,v in enumerate([item for sublist in partner2_bfactor for item in sublist]) }

    p1_df = pd.DataFrame.from_dict(p1_dict, orient='index', columns=['residue', 'pLDDT'])
    p2_df = pd.DataFrame.from_dict(p2_dict, orient='index', columns=['residue', 'pLDDT'])

    if output:
        p1_df.to_csv(outfile+'_p1.csv', index=False)
        p2_df.to_csv(outfile+'_p2.csv', index=False)
        
    return p1_df, p2_df
