
import autotst
from autotst.conformer.brute_force import *
from autotst.conformer.ga import *
from autotst.conformer.simple_es import *
from autotst.conformer.utilities import *
from hotbit.aseinterface import Hotbit
from time import time

mols = [
    "CC",
    "CCC",
    "CCCC",
    "CCCCC",
    "CCCCCC",
    "CCCCCCC",
    "CCCCCCCC",
    "CCCCCCCCC",
    "CCCCCCCCCC",
    "CCCCCCCCCCC"
]

total_trials = len(mols) * 3 * 10 # 10 of each of the 3 methods

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
smiles = mols[i % len(mols)]
ith = (i - (i % len(mols))) / len(mols)
logging.info("Job number {0} is the {1}th optimization of {2}.".format(job_number, ith, smiles))

autotst_mol = AutoTST_Molecule(smiles)

autotst_mol.ase_molecule.set_calculator(Hotbit())

if i < total_trials / 3:
    logging.info("Performing GA")
    t1 = time()
    results, conformers = perform_ga(autotst_mol, store_results=False)
    t2 = time()

elif i < (2 * total_trials / 3):
    logging.info("Performing BF")
    t1 = time() 
    results, conformers = perform_brute_force(autotst_mol, store_results=False) 
    t2 = time()

else:
    logging.info("Performing ES")
    t1 = time()
    results, conformers = perform_simple_es(autotst_mol, store_results=False)
    t2 = time()

logging.info("We have identified the following conformers from the above method:")
print conformers

logging.info("This result took: {} s to run".format(t2 - t1))

print ith, smiles, t2 - t1
