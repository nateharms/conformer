time_results = []
from time import time
from autotst.molecule import *
from autotst.conformer.brute_force import *
from ase.calculators.emt import EMT
from autotst.conformer.ga import *
from autotst.conformer.simple_es import *
import pandas as pd

smiles = [
    "CCCC", #1
    "CCCCC", #2
    "CCCCCC", #3
    "CCCCCCC", #4
    "CCCCCCCC", #5
    #"CCCCCCCCC", #6
    #"CCCCCCCCCC", #7
    #"CCCCCCCCCCC" #8
]

for i in range(5):
    for s in smiles:
        mol = AutoTST_Molecule(s)
        mol.ase_molecule.set_calculator(EMT())

        t_0 = time()

        perform_brute_force(mol, store_results=False)
        delta_bf = time() - t_0

        t_0 = time()
        perform_ga(mol)
        delta_ga = time() - t_0

        t_0 = time()
        perform_simple_es(mol)
        delta_es = time() - t_0

        num = len(s) - 3

        time_results.append([num, delta_bf, delta_ga, delta_es])
    
df = pd.DataFrame(time_results, columns=["num", "bf", "ga", "es"])
df.to_csv("profiling.csv")

