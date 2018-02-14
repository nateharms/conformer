import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from multi_molecule import Multi_Molecule

import ase
import os
import sys
from ase import io
from ase.visualize import view
import pandas as pd
from ase.calculators.gaussian import Gaussian
from ase.optimize import BFGS
import logging
FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)
#%matplotlib inline
#import matplotlib.pyplot as plt

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

i = job_number - 1


master = os.listdir("/home/harms.n/Code/ga_conformer/master/")


m = master[i]
name = m.split("_")[0]
am = io.read("/home/harms.n/Code/ga_conformer/master/" + m)
print "Reading in the following file: {}".format(m)

label = "{}_master".format(name)

calc = Gaussian(mem="5GB", nprocshared="20", label=label, scratch="/gss_gpfs_scratch/harms.n/drug_conformer", method="m062x", basis="6-311+g(2df,2p)")

am.set_calculator(calc)
opt = BFGS(am)
opt.run()

am.write("~/Code/ga_conformer/master/{0}_master_optimized.xyz".format(name, job_number))
