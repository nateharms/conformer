import os
import sys
import logging
FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)

import itertools
import random
import numpy as np
from numpy import array
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import cPickle as pickle


# do this before we have a chance to import openbabel!
import rdkit, rdkit.Chem.rdDistGeom, rdkit.DistanceGeometry

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase

import py3Dmol

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction


from multi_molecule import *
from multi_reaction import *

from ase.calculators.morse import * #chosing this calculator for now because it's fast
from ase.calculators.dftb import *
from ase.calculators.lj import *
from ase.calculators.emt import *

from copy import deepcopy

def create_initial_population(multi_object, calc=None, delta=30, population_size=200):
    """
    A function designed to take a multi_molecule, multi_rxn or multi_ts object
    and create an initial population of conformers.

    :param:
     multi_object: a multi_molecule, multi_ts, or multi_rxn object
     calc: an ASE calculator. If none is chosen, an EMT() calculator will be used
     delta: the step size between possible dihedral angles in degrees.
     population_size: the number of individuals to be used for the population

    :return:
     df: a DataFrame of the results sorted by the lowest energy conformers
    """

    df = None

    if not calc:
        calc = EMT()

    possible_dihedrals = np.arange(0, 360, delta)
    population = []

    if "Multi_Molecule" in str(multi_object.__class__):

        mol = multi_object
        mol.ase_molecule.set_calculator(calc)

        for indivudual in range(population_size):
            dihedrals = []
            for torsion in mol.torsions:
                dihedral = np.random.choice(possible_dihedrals)
                dihedrals.append(dihedral)
                i, j, k, l = torsion.indices
                RHS = torsion.RHS

                mol.ase_molecule.set_dihedral(
                    a1 = i,
                    a2 = j,
                    a3 = k,
                    a4 = l,
                    angle = float(dihedral),
                    indices = RHS
                )

            mol.update_geometry_from_ase_mol()

            e = mol.ase_molecule.get_potential_energy()

            population.append([e] + dihedrals)

    elif "Multi_Reaction" in str(multi_object.__class__):

        TS = multi_object.multi_ts

        TS.ase_ts.set_calculator(calc)

        for i in range(population_size):
            dihedrals = []

            for torsion in TS.torsions:
                dihedral = np.random.choice(possible_dihedrals)
                dihedrals.append(dihedral)
                i, j, k, l = torsion.indices
                RHS = torsion.RHS

                TS.ase_ts.set_dihedral(
                    a1=i,
                    a2=j,
                    a3=k,
                    a4=l,
                    angle=float(dihedral),
                    indices=RHS
                )

            TS.update_ts_from_ase_ts()

            e = TS.ase_ts.get_potential_energy()

            population.append([e] + dihedrals)

    elif "Multi_TS" in str(multi_object.__class__):
        TS = multi_object

        TS.ase_ts.set_calculator(calc)

        for i in range(population_size):
            dihedrals = []

            for torsion in TS.torsions:
                dihedral = np.random.choice(possible_dihedrals)
                dihedrals.append(dihedral)
                i, j, k, l = torsion.indices
                RHS = torsion.RHS

                TS.ase_ts.set_dihedral(
                    a1=i,
                    a2=j,
                    a3=k,
                    a4=l,
                    angle=float(dihedral),
                    indices=RHS
                )

            TS.update_ts_from_ase_ts()

            e = TS.ase_ts.get_potential_energy()

            population.append([e] + dihedrals)


    if len(population) > 0:
        df = pd.DataFrame(population)
        columns = ["Energy"]
        for i in range(len(multi_object.torsion_list)):
            columns = columns + ["Torsion " + str(i)]
        df.columns = columns
        df = df.sort_values("Energy")

    return df

def select_top_population(df=None, top_percent=0.30):
    """
    :param:
     df: a DataFrame of a population of torsions with columns of `Energy` and `Torsion N`
     top_percent: a float of the top percentage of the population requested

    :return:
     top: a DataFrame containing the top percentage of the population
    """

    population_size = df.shape[1]
    top_population = population_size * top_percent

    top = df.iloc[:int(top_population), :]

    return top

def perform_ga(multi_object,
               df,
               top_percent=0.3,
               tolerance = 1e-4,
               max_generations=0,
               store_generations=False,
               mutation_probability=0.2,
               delta=30):

    """
    Performs a genetic algorithm to determine the lowest energy conformer of a TS or molecule

    :param multi_object: a multi_ts, multi_rxn, or multi_molecule that you want to perform conformer analysis on
       * the ase_object of the multi_object must have a calculator attached to it.
    :param df: a DataFrame containing the initial population
    :param top_percent: float of the top percentage of conformers you want to select
    :param tolerance: float of one of the possible cut off points for the analysis
    :param max_generations: int of one of the possible cut off points for the analysis
    :param store_generations: do you want to store pickle files of each generation
    :param mutation_probability: float of the chance of mutation
    :param delta: the degree change in dihedral angle between each possible dihedral angle
    :return df: a DataFrame containing the final generation
    """
    possible_dihedrals = np.arange(0, 360, delta)
    top = select_top_population(df,
                                top_percent=top_percent
                                )

    top_population = top.shape[1]

    parent_0, parent_1 = random.sample(np.arange(top.shape[1]), 2)

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

                if random.random() < mutation_probability:
                    dihedral = np.random.choice(possible_dihedrals)
                else:
                    if 0.5 > random.random():
                        dihedral = df.iloc[parent_0, index + 1]

                    else:
                        dihedral = df.iloc[parent_1, index + 1]



                i, j, k, l = torsion.indices
                RHS = torsion.RHS

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
            generation_name = "ga_generation_{}.pkl".format(generation_tracker)

            f = open(generation_name, "w")

            pickle.dump(df, f)
        top = df.iloc[:int(top_population), :]

    return df





