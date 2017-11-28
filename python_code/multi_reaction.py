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

import rdkit, rdkit.Chem.rdDistGeom, rdkit.DistanceGeometry

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase
from rdkit.Chem.rdMolTransforms import *
from rdkit.Chem.rdChemReactions import ChemicalReaction

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

from multi_molecule import *


rmg_database = RMGDatabase()
database_path = os.path.abspath(os.path.join(os.getenv('RMGpy', '..'), '..', 'RMG-database', 'input'))
rmg_database.load(database_path,
                 kineticsFamilies=['H_Abstraction'],
                 transportLibraries=[],
                 reactionLibraries=[],
                 seedMechanisms=[],
                 thermoLibraries=['primaryThermoLibrary', 'KlippensteinH2O2', 'thermo_DFT_CCSDTF12_BAC', 'CBS_QB3_1dHR' ],
                 solvation=False,
                 )

f = open("../ts_database.pkl", "r")
ts_database = pkl.load(f)

settings = QMSettings(
    software='gaussian',
    method='m062x',
    fileStore=os.path.expandvars('.'),
    scratchDirectory=os.path.expandvars('.'),
    )

class Multi_Reaction():

    def __init__(self, label, reaction_family):
        self.label = label
        self.reaction_family = reaction_family
        self.get_reactants_and_products()
        self.get_rmg_reaction()


    def get_reactants_and_products(self):
        reactants, products = self.label.split("_")

        if "+" in reactants:
            reactants = reactants.split("+")

        if "+" in products:
            products = products.split("+")

        reactant_mols = []
        product_mols = []

        for reactant in reactants:
            reactant_mols.append(Multi_Molecule(reactant))

        for product in products:
            product_mols.append(Multi_Molecule(product))

        self.reactant_mols = reactant_mols
        self.product_mols = product_mols

    def get_rmg_reaction(self):
        """
        This method creates a labeled rmg_reaction from the reaction string
        """

        rmg_reactants = []
        rmg_products = []

        for reactant_mol in self.reactant_mols:
            rmg_reactants.append(reactant_mol.rmg_molecule)

        for product_mol in self.product_mols:
            rmg_products.append(product_mol.rmg_molecule)

        test_reaction = Reaction(reactants=rmg_reactants, products=rmg_products, reversible=True)

        reaction_list = rmg_database.kinetics.generateReactionsFromFamilies(
            reactants,
            products,
            only_families=self.reaction_family)

        for reaction in reaction_list:
            # Check if any of the RMG proposed reactions matches the reaction in the mechanism
            if test_reaction.isIsomorphic(reaction):
                atom_labels_reactants = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])
                atom_labels_products = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])

                for reactant in reaction.reactants:
                    reactant.clearLabeledAtoms()
                    for atom in reactant.atoms:
                        for atom_label in reaction.labeledAtoms:
                            if atom == atom_label[1]:
                                atom.label = atom_label[0]
                                atom_labels_reactants[atom_label[0]] = True

                for product in reaction.products:
                    product.clearLabeledAtoms()
                    for atom in product.atoms:
                        for atom_label in reaction.labeledAtoms:
                            if atom == atom_label[1]:
                                atom.label = atom_label[0]
                                atom_labels_products[atom_label[0]] = True

                if all(atom_labels_reactants.values()) and all(atom_labels_products.values()):
                    # We successfully labeled all of the atoms
                    break
        self.rmg_reaction = reaction
        self.rmg_qm_reaction = QMReaction(reaction=reaction, settings=settings, tsDatabase=ts_database)





Multi_Reaction(label)