import os
import sys
import cPickle as pkl
import logging
FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)

import re
import imp
import itertools
import random
import numpy as np
from numpy import array
import pandas as pd

# do this before we have a chance to import openbabel!
import rdkit, rdkit.Chem.rdDistGeom, rdkit.DistanceGeometry

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit import rdBase
from rdkit.Chem.rdMolTransforms import *

import py3Dmol

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.kinetics import PDepArrhenius, PDepKineticsModel

from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.kinetics import KineticsDepository, KineticsRules
from rmgpy.qm.main import QMCalculator, QMSettings
from rmgpy.qm.qmdata import QMData
from rmgpy.qm.reaction import QMReaction
from rmgpy.qm.molecule import QMMolecule

rxnFamily = ['H_Abstraction'] # Can be set to other reaction families as well
print('Loading RMG Database ...')
rmgDatabase = RMGDatabase()
databasePath = os.path.abspath(os.path.join(os.getenv('RMGpy', '..'), '..', 'RMG-database', 'input'))
print(databasePath)
rmgDatabase.load(databasePath,
                 kineticsFamilies=['H_Abstraction'],
                 transportLibraries=[],
                 reactionLibraries=[],
                 seedMechanisms=[],
                 thermoLibraries=['primaryThermoLibrary', 'KlippensteinH2O2', 'thermo_DFT_CCSDTF12_BAC', 'CBS_QB3_1dHR' ],
                 solvation=False,
                 )
print('RMG Database Loaded')


## We are using a ts Pickle from discovery.
## I was running into issues performing ts opts locally, so this is a workaround
f = open("../tsDatabase.pkl", "r")
tsDatabase = pkl.load(f)

settings = QMSettings(
    software='gaussian',
    method='m062x',
    fileStore=os.path.expandvars('.'),
    scratchDirectory=os.path.expandvars('.'),
    )

def get_ts_rdmol(input_reaction):
    """

    """
    rxnFamilies = ["H_Abstraction"]
    rSpecies1, rSpecies2 = input_reaction.reactants
    pSpecies1, pSpecies2 = input_reaction.products

    for species in itertools.chain(input_reaction.reactants, input_reaction.products):
            species = species.generateResonanceIsomers()

    testReaction = Reaction(
        reactants = input_reaction.reactants,
        products = input_reaction.products,
        reversible = True)

    reactants = [species for species in input_reaction.reactants]
    products = [species for species in input_reaction.products]

    reactionList = []

    checkRxn = rmgDatabase.kinetics.generateReactionsFromFamilies(
        reactants,
        products,
        only_families=rxnFamilies)

    for rxn0 in checkRxn:
        reactionList.append(rxn0)
    logging.info("generateReactionsFromFamilies successfuly performed")

    gotOne=False
    for reaction in reactionList:
        # Check if any of the RMG proposed reactions matches the reaction in the mechanism
        if testReaction.isIsomorphic(reaction):
            #print "Found matching reaction"
            # Now add the labeled atoms to the Molecule, and check all labels were added
            atLblsR = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])
            atLblsP = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])

            for reactant in reaction.reactants:
                reactant.clearLabeledAtoms()
                for atom in reactant.atoms:
                    for atomLabel in reaction.labeledAtoms:
                        if atom==atomLabel[1]:
                            atom.label = atomLabel[0]
                            atLblsR[atomLabel[0]] = True

            for product in reaction.products:
                product.clearLabeledAtoms()
                for atom in product.atoms:
                    for atomLabel in reaction.labeledAtoms:
                        if atom==atomLabel[1]:
                            atom.label = atomLabel[0]
                            atLblsP[atomLabel[0]] = True

            if all( atLblsR.values() ) and all( atLblsP.values() ):
                # We successfully labeled all of the atoms
                gotOne=True
                break

    rxn = QMReaction(reaction=reaction, settings=settings, tsDatabase=tsDatabase)


    reactant, product = rxn.setupMolecules()
    mol = reactant.toRDKitMol(removeHs=False)

    AllChem.EmbedMolecule(mol)
    labels, atomMatch = rxn.getLabels(reactant)
    tsRDMol, bm, rxn.reactantGeom = rxn.generateBoundsMatrix(reactant)
    bm = rxn.editMatrix(reactant, bm, labels)
    tsRDMol = rxn.reactantGeom.rd_embed(tsRDMol, 15, bm=bm, match=atomMatch)[0]
    tsRDMol = rdkit.Chem.rdchem.RWMol(tsRDMol)

    return tsRDMol

def create_pseudo_geometry(tsRDMol):
    """
    Creates a `pseduo` bond between the reaction center molecule.
    This allows for different torsions to be applied about the reaction center as well.

    Input:
    * tsRDMol: rdkit.mol object, a TS geometry of two seperate molecules

    Output:
    * tsRDMol: rdkit.mol object, a TS geometry of a `single` molecule
    """

    for atom in tsRDMol.GetAtoms():
        idx = atom.GetIdx()
        num = atom.GetAtomicNum()
        rmg_atom = reactant.atoms[idx]

        if rmg_atom.label:
            if rmg_atom.label == "*1":
                atom1_star = atom
            if rmg_atom.label == "*2":
                atom2_star = atom
            if rmg_atom.label == "*3":
                atom3_star = atom

    bond_between_23 = False
    try:
        tsRDMol.AddBond(atom1_star.GetIdx(), atom2_star.GetIdx(), order=rdkit.Chem.rdchem.BondType.SINGLE)
    except RuntimeError:
        #print "Bond already exists betwee 1* and 2*"
        bond_between_23 = True
        tsRDMol.AddBond(atom2_star.GetIdx(), atom3_star.GetIdx(), order=rdkit.Chem.rdchem.BondType.SINGLE)

    return tsRDMol


def get_ts_torsion_energies(tsRDMol, delta):

    torsion_list = get_torsion_list(tsRDMol)

    torsion_combos = list( itertools.combinations_with_replacement(torsion_angles, len(torsion_list)) )
    if len(torsion_list) != 1:
        torsion_combos = list(
            set(
                torsion_combos +
                list(itertools.combinations_with_replacement(
                    torsion_angles[::-1], len(torsion_list)
                ))))
    df = []

    tup = tsRDMol.GetConformers()
    conformer = tup[0]

    for combo in torsion_combos:
        geometry = zip(torsion_list, combo)
        for torsion in geometry:
            #print torsion
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

        if bond_between_23 == False:
            tsRDMol.RemoveBond(atom1_star.GetIdx(), atom2_star.GetIdx())
        else:
            tsRDMol.RemoveBond(atom2_star.GetIdx(), atom3_star.GetIdx())

        e0, m0 = calc_energy_ts(tsRDMol)
        e1, m1 = constrained_energy_calc_ts(tsRDMol, torsion_list, combo)
        e2, m2 = optimize_and_calc_ts(tsRDMol, torsion_list, combo)

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


def get_lowest_ts_energy_conformers(tsRDMol, delta):
    """
    A function that returns a dict of the lowest energy conformers calculated in a variety of different ways

    Inputs:
    * tsRDMol: a rdkit molecule object of your TS
    * delta: `int` that describes the increment for the scan

    Outputs:
    * conformer_dict: `dict` containing optimization energies and molecules corresponding to
                     three differnet calculators

    """

    conformer_dict = {}

    df = get_ts_torsion_energies(tsRDMol, delta)

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
