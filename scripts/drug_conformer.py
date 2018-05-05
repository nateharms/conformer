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

import rdkit, rdkit.Chem.rdDistGeom, rdkit.DistanceGeometry
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction

from ase.optimize import BFGS
from ase.calculators.emt import *
from ase.calculators.gaussian import *

from autotst.molecule import *
from autotst.conformer.utilities  import *
from autotst.conformer.ga import *
from autotst.conformer.simple_es import *
from autotst.conformer.brute_force import *

if len(sys.argv)>1:
    job_number = int(sys.argv[-1])
elif os.getenv('SLURM_ARRAY_TASK_ID'):
    job_number = int(os.getenv('SLURM_ARRAY_TASK_ID'))
elif os.getenv('LSB_JOBINDEX'):
    job_number = int(os.getenv('LSB_JOBINDEX'))
else:
    #raise Exception("Specif y a TS number!")
    logging.warning("Number not specified as script argument or via environment variable, so using default")
    job_number = 1

logging.info("RUNNING WITH JOB NUMBER i = {}".format(job_number))

# For this test, we will be looking at Paclitaxel. A commonly used anti-cancer drug
mol = AutoTST_Molecule("CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C")

# Setting the population size to 25 times the numer of torsions. This just seems reasonable
population_size = len(mol.torsions) * 25

# For the initial calculations, we are going to use EMT to identify the lowest energy conformer,
# then we will optimize with gaussian.

mol.ase_molecule.set_calculator(EMT())

initial_pop = create_initial_population(multi_object=mol, delta=float(30), population_size=population_size)

if job_number == 1:
    results = perform_brute_force(mol,
                            delta=float(30),
                            store_directory=".")

    results.to_csv("brute_force_results.csv")



    for h, tor in enumerate(brute_force_mol.torsions):
        dihedral = results.iloc[0][h+1]
        i,j,k,l = tor.indices
        RHS = tor.right_mask
        mol.ase_molecule.set_dihedral(a1 = i,
                                      a2 = j,
                                      a3 = k,
                                      a4 = l,
                                      angle= float( dihedral),
                                      mask=RHS)


    mol.ase_molecule.set_calculator(Gaussian(label="BF"))
    opt = BFGS(atoms=mol.ase_molecule)
    opt.run(fmax=0.05)

    mol.update_from_ase_mol()

    print "The optomized geometry from brute_force is:"

    for atom in mol.rmg_molecule.atoms:
        print("{0} \t {1}".format(atom, atom.coords))

    Chem.MolToMolFile(mol.rdkit_molecule, 'brute_force_paclitaxel.mol')

if job_number == 2:

    results = perform_ga(mol,
               initial_pop,
               top_percent=0.3,
               tolerance=1e-4,
               store_generations=True,
               store_directory=".",
               mutation_probability=0.2,
               delta=float(30))

    for h, tor in enumerate(mol.torsions):
        dihedral = results.iloc[0][h+1]
        i,j,k,l = tor.indices
        RHS = tor.right_mask
        mol.ase_molecule.set_dihedral(a1 = i,
                                         a2 = j,
                                         a3 = k,
                                         a4 = l,
                                         angle= float(dihedral),
                                         mask=RHS)
    mol.ase_molecule.set_calculator(Gaussian(label="GA"))
    opt = BFGS(atoms=mol.ase_molecule)
    opt.run(fmax=0.05)

    print "The optomized geometry from GA is:"

    for atom in mol.rmg_molecule.atoms:
        print("{0} \t {1}".format(atom, atom.coords))

    Chem.MolToMolFile(mol.rdkit_molecule, 'ga_paclitaxel.mol')

if job_number == 3:

    results = perform_simple_es(mol,
                      initial_pop,
                      top_percent=0.3,
                      tolerance=1e-5,
                      max_generations=100,
                      store_generations=True,
                      store_directory=".")


    for h, tor in enumerate(mol.torsions):
        dihedral = results.iloc[0][h+1]
        i,j,k,l = tor.indices
        RHS = tor.right_mask
        mol.ase_molecule.set_dihedral(a1 = i,
                                         a2 = j,
                                         a3 = k,
                                         a4 = l,
                                         angle= float(dihedral),
                                         mask=RHS)


    mol.ase_molecule.set_calculator(Gaussian(label="ES"))
    opt = BFGS(atoms=mol.ase_molecule)
    opt.run(fmax=0.05)
    print "The optomized geometry from ES is:"

    for atom in mol.rmg_molecule.atoms:
        print("{0} \t {1}".format(atom, atom.coords))

    Chem.MolToMolFile(mol.rdkit_molecule, 'es_paclitaxel.mol')
