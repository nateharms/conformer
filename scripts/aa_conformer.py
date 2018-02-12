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

from ase.calculators.emt import *
from ase.calculators.gaussian import *
from simple_es import *
from utilities import *

from ase.optimize import BFGS


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
#i = i + 999 #### ADDED THIS LINE TO GET ARRAY OVER 1000 FOR SLURM
logging.info("RUNNING WITH JOB NUMBER i = {}".format(job_number))


smiles_dict = {
    "Phe" : "NC(CC1=CC=CC=C1)C(=O)O",
    "Gly" : "NCCC(=O)O",
    "Ala" : "NC(C)C(=O)O",
    "Ile" : "NC(C(C)CC)C(=O)O",
    "Leu" : "NC(CC(C)C)C(=O)O",
    "Val" : "NC(C(C)C)C(=O)O",
    "Trp" : "NC(CC1=CNC2=C1(C=CC=C2))C(=O)O"
}

i = job_number - 1
name, smiles = list(smiles_dict.iteritems())[i % 7]
print "Job number {0} is the optimization of {1}.".format(job_number, name)

if not os.path.exists("/home/harms.n/Code/ga_conformer/results/{0}_{1}.xyz".format(name, job_number)):
    print "New job, running it"

    mol = Multi_Molecule(smiles)
    mol.ase_molecule.set_calculator(EMT())
    initial_pop = create_initial_population(mol)
    top_pop = perform_simple_es(mol, initial_pop)

    for i, dihedral in enumerate(top_pop.iloc[0,1:]):
        tor = mol.torsions[i]
        mol.ase_molecule.set_dihedral(tor.indices, angle=dihedral, mask=tor.right_mask)
    mol.update_geometry_from_ase_mol()
    label = "{0}_{1}".format(name, job_number)#str(name) + "_" + str(job_number)
    calc = Gaussian(mem="5GB", nprocshared="20", label=label, scratch="/gss_gpfs_scratch/harms.n/drug_conformer", method="m062x", basis="6-311+g(2df,2p)")
    mol.ase_molecule.set_calculator(calc = calc)
    opt = BFGS(mol.ase_molecule)
    opt.run()

    mol.ase_molecule.write("~/Code/ga_conformer/results/{0}_{1}.xyz".format(name, job_number))

else:
    print "Previous job complete. See the below file."
    print "~/Code/ga_conformer/results/{0}_{1}.xyz".format(name, job_number)
print "Job complete!!!"
