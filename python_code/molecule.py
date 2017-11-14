import os
import sys

import re
import imp
import itertools
import random
import numpy as np
from numpy import array
import pandas as pd

import rdkit, rdkit.Chem.rdDistGeom, rdkit.DistanceGeometry

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit import rdBase
from rdkit.Chem.rdMolTransforms import *

from calculator import *



def get_rdmol(SMILES):
    """
    A function to get a rewritable RDKit molecule object from a smiles

    This will only work for actual molecules

    Input:
    * SMILES string

    Output:
    * Embeded rewritable RDKit molecule object with hydrogens

    """
    RDMol = Chem.AddHs(Chem.MolFromSmiles(
            SMILES
        ))
    rdkit.Chem.AllChem.EmbedMolecule(RDMol)

    RDMol = rdkit.Chem.rdchem.RWMol(RDMol)

    return RDMol

def get_torsion_list(RDMol):
    """
    A function to return a list of the possible torsions in an RDKit molecule

    Input:
    * RDKit molecule

    Output:
    * list of tuples containing 4 indicies of atoms representing a torsion
    """
    torsion_list = []
    for bond1 in RDMol.GetBonds():
        atom1 = bond1.GetBeginAtom()
        atom2 = bond1.GetEndAtom()
        if atom1.IsInRing() or atom2.IsInRing():
            # Making sure that bond1 we're looking at are in a ring
            continue

        bond_list1 = list(atom1.GetBonds())
        bond_list2 = list(atom2.GetBonds())

        if not len(bond_list1) > 1 and not len(bond_list2) > 1:
            # Making sure that there are more than one bond attached to
            # the atoms we're looking at
            continue

        # Getting the 0th and 3rd atom and insuring that atoms
        # attached to the 1st and 2nd atom are not terminal hydrogens
        # We also make sure that all of the atoms are properly bound together

        # If the above are satisified, we append a tuple of the torsion our torsion_list
        got_atom0 = False
        got_atom3 = False

        for bond0 in bond_list1:
            atomX = bond0.GetOtherAtom(atom1)
            if atomX.GetAtomicNum() == 1 and len(atomX.GetBonds()) == 1:
                # This means that we have a terminal hydrogen, skip this
                # NOTE: for H_abstraction TSs, a non teminal H should exist
                continue
            if atomX.GetIdx() != atom2.GetIdx():
                got_atom0 = True
                atom0 = atomX

        for bond2 in bond_list2:
            atomY = bond2.GetOtherAtom(atom2)
            if atomY.GetAtomicNum() == 1 and len(atomY.GetBonds()) == 1:
                # This means that we have a terminal hydrogen, skip this
                continue
            if atomY.GetIdx() != atom1.GetIdx():
                got_atom3 = True
                atom3 = atomY

        if not (got_atom0 and got_atom3):
            # Making sure atom0 and atom3 were not found
            continue

        # Looking to make sure that all of the atoms are properly bonded to eached
        if (
            RDMol.GetBondBetweenAtoms(atom0.GetIdx(), atom1.GetIdx()) and
            RDMol.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx()) and
            RDMol.GetBondBetweenAtoms(atom2.GetIdx(), atom3.GetIdx())   ) :

            torsion_tup = (atom0.GetIdx(), atom1.GetIdx(), atom2.GetIdx(), atom3.GetIdx())
            torsion_list.append(torsion_tup)

    return torsion_list

def get_torsion_energies(RDMol, delta):
    """
    A function to take a rdkit molecule, identify the torsions,
    and calculate the constrained, partially optimized and fully optimized energies
    at a number of torsions ranging from 0 to 360 degrees at `delta` degree increments

    Input:
    * RDMol: A rdkit molecule object
    * delta: `int` that describes the increment for the scan

    Output:
    * df: a pandas dataframe containing the following columns
        * Constrained energy
        * The rdkit molecule for which the constrained energy was calculated for
        * Partially optimized energy
        * The rdkit molecule for which the partially optimized energy was calculated for
        * Fully optimized energy
         The rdkit molecule for which the fully optimized energy was calculated for
        * Torsional angles for each possible torsion (starting at 0)

    """

    torsion_list = get_torsion_list(RDMol)

    torsion_angles = np.arange(0, 360 + delta, delta)
    torsion_combos = list( itertools.combinations_with_replacement(torsion_angles, len(torsion_list)) )
    if len(torsion_list) != 1:
        torsion_combos = list(
            set(
                torsion_combos +
                list(itertools.combinations_with_replacement(
                    torsion_angles[::-1], len(torsion_list)
                ))))

    df = []

    tup = RDMol.GetConformers()
    conformer = tup[0]
    for combo in torsion_combos:


        geometry = zip(torsion_list, combo)

        for torsion in geometry:

            i = torsion[0][0]
            j = torsion[0][1]
            k = torsion[0][2]
            l = torsion[0][3]
            angle = torsion[1]

            SetDihedralDeg(conformer,
                           i,
                           j,
                           k,
                           l,
                           angle)

        e0, m0 = calc_energy(RDMol)
        e1, m1 = constrained_energy_calc(RDMol, torsion_list, combo)
        e2, m2 = optimize_and_calc(RDMol, torsion_list, combo)

        df.append( [e0, e1, e2,
                   m0, m1, m2]
                  + list(combo))


    df = pd.DataFrame(df)
    columns = ["No Optimization Energy", "Partial Optimization Energy", "Full Optimization Energy",
              "No Optimization Mol", "Partial Optimization Mol", "Full Optimization Mol"]
    for i in range(len(torsion_list)):
        columns.append("Torsion " + str(i))

    df.columns = columns

    return df

def get_lowest_energy_conformers(RDMol, delta):
    """
    A function that returns a dict of the lowest energy conformers calculated in a variety of different ways

    Inputs:
    * RDMol: a rdkit molecule object
    * delta: `int` that describes the increment for the scan

    Outputs:
    * conformer_dict: `dict` containing optimization energies and molecules corresponding to
                     three differnet calculators

    """

    conformer_dict = {}

    df = get_torsion_energies(RDMol, delta)

    conformer_dict["No Optimization"] = (
        df.sort("No Optimization Energy")["No Optimization Energy"].iloc[0],
        df.sort("No Optimization Energy")["No Optimization Mol"].iloc[0]
    )

    conformer_dict["Partial Optimization"] = (
        df.sort("Partial Optimization Energy")["Partial Optimization Energy"].iloc[0],
        df.sort("Partial Optimization Energy")["Partial Optimization Mol"].iloc[0]
    )

    conformer_dict["Full Optimization"] = (
        df.sort("Full Optimization Energy")["Full Optimization Energy"].iloc[0],
        df.sort("Full Optimization Energy")["Full Optimization Mol"].iloc[0]
    )

    return conformer_dict
