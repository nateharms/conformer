import os
import sys
import ase
from ase.io import read
from ase.io.xyz import write_xyz
from ase.calculators.gaussian import Gaussian
from ase.optimize import BFGS
import logging
import numpy as np

FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)


index = int(sys.argv[-1])

if os.getenv('SLURM_ARRAY_TASK_ID'):
    job_number = int(os.getenv('SLURM_ARRAY_TASK_ID'))
elif os.getenv('LSB_JOBINDEX'):
    job_number = int(os.getenv('LSB_JOBINDEX'))
else:
    #raise Exception("Specif y a TS number!")
    logging.warning("Number not specified as script argument or via environment variable, so using default")
    job_number = 1

logging.info("RUNNING WITH JOB NUMBER i = {}".format(job_number))

i = job_number - 1

path = "/gss_gpfs_scratch/harms.n/conformers/new_ga/{}"

mols_of_interest = [
    "Ala",
    "Gly",
    "Ile",
    "Leu",
    "Phe",
    "Trp",
    "Val"
]

if index == 8:
    mol = "Val"
    i += 1000
elif index > 8:
    mol = "Leu"
    i += (index - 8) * 1000
else:
    mol = mols_of_interest[index - 1]

print "This is the {}th optimization of {}...".format(i, mol)
files_to_optimize = []  
list_dir =  os.listdir(path.format(mol))
for f in list_dir:
    if f.endswith(".xyz"):
        files_to_optimize.append(os.path.join(path.format(mol), f))

files_to_optimize = np.array(files_to_optimize)
files_to_optimize.sort()
file_of_interest = files_to_optimize[i]

print "We are optimizing this file:\t{}".format(file_of_interest)

atoms = read(file_of_interest)

scratch = "/gss_gpfs_scratch/harms.n/conformers/"

ff = file_of_interest.split("/")[-1]
label = ff.split(".")[0]


os.chdir(scratch)

calc = Gaussian(mem="5GB",
                nprocshared=20,
                label=label,
                scratch=".",
                method="m062x",
                basis="6-311+g(2df,2p)",
                multiplicity=1
                )

atoms.set_calculator(calc)
opt = BFGS(atoms=atoms)
opt.run()

results_path = os.path.join("/gss_gpfs_scratch/harms.n/conformers/new_ga_results", mol, file_of_interest)


f = open(results_path, "w")
write_xyz(f, atoms)

print 
print "The potential energy of {} is:".format(file_of_interest[20:-4].replace("/","_"))
print atoms.get_potential_energy()
