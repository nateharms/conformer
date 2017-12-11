from multi_molecule import *
from multi_reaction import *
from brute_force import *
from ga import *
from simple_es import *
from ase.calculators.emt import *
from ase.calculators.lj import *

multi_object = Multi_Reaction("OCCCC+[O]O_OCC[CH]C+OO", "H_Abstraction")
multi_object.multi_ts.ase_ts.set_calculator(EMT())