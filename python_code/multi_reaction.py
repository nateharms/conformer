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

    def __init__(self, label, reaction_family, rmg_reaction = None):
        self.label = label
        self.reaction_family = reaction_family

        if rmg_reaction:

            reactant_mols = []
            product_mols = []

            reactants = rmg_reaction.reactants
            products = rmg_reaction.products

            for reactant in reactants:
                reactant_mols.append(Multi_Molecule(reactant.SMILES))

            for product in products:
                product_mols.append(Multi_Molecule(product))

            self.reactant_mols = reactant_mols
            self.product_mols = product_mols
        else:
            self.get_reactants_and_products()

        self.get_rmg_reactions()
        self.create_ts_geometries()

    def get_reactants_and_products(self):

        """
        This uses the reaction label to creat multi_molecule objects of

        :return:
        """
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


    def get_rmg_reactions(self):
        """
        This method creates a labeled rmg_reaction from the reaction string

        the reaction string should look as follows: r1+r2_p1+p2
        """

        rmg_reactants = []
        rmg_products = []

        for reactant_mol in self.reactant_mols:
            rmg_reactants.append(reactant_mol.rmg_molecule)

        for product_mol in self.product_mols:
            rmg_products.append(product_mol.rmg_molecule)

        test_reaction = Reaction(reactants=rmg_reactants, products=rmg_products, reversible=True)

        reaction_list = rmg_database.kinetics.generateReactionsFromFamilies(
            rmg_reactants,
            rmg_products,
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



    def create_ts_geometries(self):
        """
        A method to use the tools in rmg / autotst to create a reasonable TS geometry
        This will create the geometry in both rdkit and ase

        :return:
        self.multi_ts: a multi_ts object that contains geometries of a ts in
                        rdkit, ase, and rmg molecules
        """

        self.multi_ts = Multi_TS(self)


class Multi_TS():
    def __init__(self, Multi_Reaction):

        self.multi_reaction = Multi_Reaction
        self.label = Multi_Reaction.label #make sure that both the reaction and TS have same label

        self.create_rdkit_ts_geometry()
        self.create_ase_ts_geometry()
        self.create_rmg_ts_geometry()

    def create_rdkit_ts_geometry(self):

        self.rmg_ts, _ = self.multi_reaction.rmg_qm_reaction.setupMolecules()

        labels, atom_match = self.multi_reaction.rmg_qm_reaction.getLabels(self.rmg_ts)

        rdkit_ts, bm, self.multi_reaction.rmg_qm_reaction.reactantGeom = self.multi_reaction.rmg_qm_reaction.generateBoundsMatrix(self.rmg_ts)

        bm = self.multi_reaction.rmg_qm_reaction.editMatrix(self.rmg_ts, bm, labels)

        self.rdkit_ts = self.multi_reaction.rmg_qm_reaction.reactantGeom.rd_embed(rdkit_ts, 15, bm=bm, match=atom_match)[0]

        #os.remove("*.mol") #removing unnecessary .mol files

    def create_ase_ts_geometry(self):

        mol_list = AllChem.MolToMolBlock(self.rdkit_ts).split('\n')
        ase_atoms = []
        for i, line in enumerate(mol_list):

            if i > 3:

                try:
                    atom0, atom1, bond, rest = line
                    atom0 = int(atom0)
                    atom0 = int(atom1)
                    bond = float(bond)

                except ValueError:
                    try:
                        x, y, z, symbol = line.split()[0:4]
                        x = float(x)
                        y = float(y)
                        z = float(z)
                        # print symbol

                        ase_atoms.append(Atom(symbol=symbol, position=(x, y, z)))

                    except:
                        continue

        self.ase_ts = Atoms(ase_atoms)

    def create_rmg_ts_geometry(self):

        for i, position in enumerate(self.ase_ts.get_positions()):
            self.rmg_ts.atoms[i].coords = position

