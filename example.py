# A test file to act as a tutorial of how this tool works

from multi_molecule import *
from multi_reaction import *
from brute_force import *
from ga import *
from simple_es import *
from ase.calculators.emt import *

# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Performing GA and Simple-ES on a molecule

mol = Multi_Molecule("CCCOC")
mol.ase_molecule.set_calculator(EMT())
initial_pop = create_initial_population(multi_object=mol)

"""perform_brute_force(mol,
                    delta=float(30),
                    store_directory="./example_results")

perform_ga(mol,
           initial_pop,
           top_percent=0.3,
           tolerance=1e-5,
           max_generations=5,
           store_generations=True,
           store_directory="./example_results",
           mutation_probability=0.2,
           delta=30)

perform_simple_es(mol,
                  initial_pop,
                  top_percent=0.3,
                  tolerance=1e-5,
                  max_generations=5,
                  store_generations=True,
                  store_directory="./example_results")"""

# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Performing GA and Simple-ES on a reaction

rxn = Multi_Reaction("OCCCC+[O]O_OCC[CH]C+OO", "H_Abstraction")
rxn.multi_ts.ase_ts.set_calculator(EMT())
rxn_initial_pop = create_initial_population(multi_object=rxn)

"""perform_brute_force(rxn,
                    delta=float(30),
                    store_directory="./example_results")"""

"""perform_ga(rxn,
           rxn_initial_pop,
           top_percent=0.3,
           tolerance = 1e-4,
           max_generations=5,
           store_generations=True,
           store_directory="./example_results",
           mutation_probability=0.2,
           delta=30)"""

perform_simple_es(rxn,
                  rxn_initial_pop,
                  top_percent=0.3,
                  tolerance=1e-4,
                  max_generations=5,
                  store_generations=True,
                  store_directory=".")
