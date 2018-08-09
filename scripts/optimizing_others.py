import os
import sys
from autotst.molecule import *
import ase
from ase.io import read
from ase.calculators.gaussian import Gaussian
from ase.calculators.nwchem import NWChem
from ase.optimize import BFGS
import logging

FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)


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

i = job_number - 1

path = "../external_databae/{}/uncapped/bare"

mols_of_interest = [
    "ala",
    "gly",
    "ile",
    "leu",
    "phe",
    "trp",
    "val"
]

files_to_optimize = []
for mol in mols_of_interest:
    
    list_dir =  os.listdir(path.format(mol))
    for f in list_dir:
        if f.endswith(".xyz"):
            files_to_optimize.append(os.path.join(path.format(mol), f))

file_of_interest = files_to_optimize[i]

atoms = read(file_of_interest)

scratch = "/projects/CPOX/northeastern_comocheng/conformers/optimizing_others/"

os.chdir(scratch)

calc = NWChem(
                label=file_of_interest[20:-4].replace("/","_"),
                scratch=scratch,
                method="m06-2x",
                basis="6-311++G(2D,2P)",
                multiplicity=1
                )

atoms.set_calculator(calc)
opt = BFGS(atoms=atoms)
opt.run()
print 
print "The potential energy of {} is:".format(file_of_interest[20:-4].replace("/","_"))
print atoms.get_potential_energy()
