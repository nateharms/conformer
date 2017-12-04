import os
import sys
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
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns


# do this before we have a chance to import openbabel!
import rdkit, rdkit.Chem.rdDistGeom, rdkit.DistanceGeometry

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase

import py3Dmol

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction

import ase
from ase import Atom, Atoms

from multi_molecule import *
from multi_reaction import *

from ase.calculators.morse import * #chosing this calculator for now because it's fast
from ase.calculators.dftb import *
from ase.calculators.lj import *
from ase.calculators.emt import *

import cPickle as pickle

from copy import deepcopy

from ga import create_initial_population, select_top_population


def perform_simple_es(multi_object,
                      df,
                      top_percent=0.3,
                      tolerance=1e-4,
                      max_generations=0,
                      store_generations=False):
    """
    Performs a simple evolutionary strategy to determine the lowest energy conformer of a TS or molecule

    :param multi_object: a multi_ts, multi_rxn, or multi_molecule that you want to perform conformer analysis on
       * the ase_object of the multi_object must have a calculator attached to it.
    :param df: a DataFrame containing the initial population
    :param top_percent: float of the top percentage of conformers you want to select
    :param tolerance: float of one of the possible cut off points for the analysis
    :param max_generations: int of one of the possible cut off points for the analysis
    :param store_generations: do you want to store pickle files of each generation
    :return df: a DataFrame containing the final generation
    """
    top = select_top_population(df,
                                top_percent=top_percent
                                )
    top_population = top.shape[1]

    population_size = df.shape[1]

    # Takes each of the molecule objects
    if "Multi_Molecule" in str(multi_object.__class__):
        ase_object = multi_object.ase_molecule
        torsions = multi_object.torsions

    elif "Multi_Reaction" in str(multi_object.__class__):
        ase_object = multi_object.multi_ts.ase_ts
        torsions = multi_object.multi_ts.torsions

    elif "Multi_TS" in str(multi_object.__class__):
        ase_object = multi_object.ase_ts
        torsions = multi_object.torsions

    gen_number = 0
    generation_tracker = 0
    while (top.std()[0] / top.mean()[0] > tolerance) or (gen_number <= max_generations):
        generation_tracker = generation_tracker + 1
        if max_generations > 0:
            # If the user sets a number of max generations, this will keep track of that and break
            # once the number of max generations is reached
            gen_number = gen_number + 1

        results = []
        for individual in range(population_size):
            dihedrals = []
            for index, torsion in enumerate(torsions):
                i, j, k, l = torsion.indices
                RHS = torsion.RHS

                dihedral = random.gauss(top.mean()[index + 1], top.std()[index + 1])
                dihedrals.append(dihedral)
                ase_object.set_dihedral(a1=i,
                                                     a2=j,
                                                     a3=k,
                                                     a4=l,
                                                     angle=float(dihedral),
                                                     indices=RHS)
                # Updating the molecule
                if "Multi_Molecule" in str(multi_object.__class__):
                    multi_object.update_geometry_from_ase_mol()

                elif "Multi_Reaction" in str(multi_object.__class__):
                    multi_object.multi_ts.update_ts_from_ase_ts()

                elif "Multi_TS" in str(multi_object.__class__):
                    multi_object.update_ts_from_ase_ts()

            e = ase_object.get_potential_energy()
            results.append([e] + dihedrals)

        df = pd.DataFrame(results)

        columns = ["Energy"]
        for i in range(len(multi_object.torsion_list)):
            columns = columns + ["Torsion " + str(i)]

        df.columns = columns
        df = df.sort_values("Energy")

        if store_generations == True:
            # This portion stores each generation if desired
            generation_name = "simple_es_generation_{}.pkl".format(generation_tracker)

            f = open(generation_name, "w")

            pickle.dump(df, f)
        top = df.iloc[:int(top_population), :]

    return df









