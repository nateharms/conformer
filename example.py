from multi_molecule import *
from multi_reaction import *
from ga import *
from simple_es import *

from ase.calculators.emt import *

mol = Multi_Molecule("CCOC")
#mol.ase_molecule.set_calculator(EMT())

initial_pop = create_initial_population(multi_object=mol)

perform_ga(mol,
           initial_pop,
           top_percent=0.3,
           #tolerance=0,
           max_generations=5,
           store_generations=True,
           mutation_probability=0.2,
           delta=30)

perform_simple_es(mol,
                  initial_pop,
                  top_percent=0.3,
                  #tolerance=0,
                  max_generations=5,
                  store_generations=True)




