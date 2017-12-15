from multi_molecule import *
from multi_reaction import *
from brute_force import *
from ga import *
from simple_es import *
from ase.calculators.emt import *
from ase.calculators.lj import *
from ase.calculators.morse import *
from ase.calculators.gaussian import *
from ase import optimize
from copy import deepcopy
import cPickle as pickle


# For this test, we will be looking at Paclitaxel. A commonly used anti-cancer drug
mol = Multi_Molecule("CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C")

# Setting the population size to 40 times the numer of torsions. This just seems reasonable
population_size = len(mol.torsions) * 40

# For the initial calculations, we are going to use EMT to identify the lowest energy conformer,
# then we will optimize with gaussian.

mol.ase_molecule.set_calculator(EMT())
initial_pop = create_initial_population(multi_object=mol, delta=float(30), population_size=population_size)


brute_force_results = perform_brute_force(multi_object,
                        delta=float(30),
                        store_directory="./drug_conformer_results"):


mol_copy = deepcopy(mol)
for h, tor in enumerate(mol_copy.torsions):
    dihedral = brute_force_results.iloc[0][h+1]
    i,j,k,l = tor.indices
    RHS = tor.right_mask
    mol_copy.ase_molecule.set_dihedral(a1 = i,
                                       a2 = j,
                                       a3 = k,
                                       a4 = l,
                                       angle= float(dihedral),
                                       mask=RHS)
mol_copy.ase_molecule.set_calculator(Gaussian())
# TODO: Need to import optimizer

opt = optimize.BFGS(atoms=mol_copy.ase_molecule)
opt.run(fmax=0.05, trajectory='./drug_conformer_results/brute_force_paclitaxel.traj')
brute_force_pickle = open("./drug_conformer_results/brute_force_final.pkl")
pickle.dump(mol_copy, brute_force_pickle)


ga_results = perform_ga(mol,
           initial_pop,
           top_percent=0.3,
           tolerance=1e-5,
           store_generations=True,
           store_directory="./drug_conformer_results/",
           mutation_probability=0.2,
           delta=float(30))

mol_copy = deepcopy(mol)
for h, tor in enumerate(mol_copy.torsions):
    dihedral = ga_results.iloc[0][h+1]
    i,j,k,l = tor.indices
    RHS = tor.right_mask
    mol_copy.ase_molecule.set_dihedral(a1 = i,
                                       a2 = j,
                                       a3 = k,
                                       a4 = l,
                                       angle= float(dihedral),
                                       mask=RHS)
mol_copy.ase_molecule.set_calculator(Gaussian())

opt = optimize.BFGS(atoms=mol_copy.ase_molecule)
opt.run(fmax=0.05, trajectory='./drug_conformer_results/ga_paclitaxel.traj')
ga_pickle = open("./drug_conformer_results/ga_final.pkl")
pickle.dump(mol_copy, ga_pickle)

simple_es_results = perform_simple_es(mol,
                  initial_pop,
                  top_percent=0.3,
                  tolerance=1e-5,
                  max_generations=100,
                  store_generations=True,
                  store_directory="./drug_conformer_results/")

mol_copy = deepcopy(mol)
for h, tor in enumerate(mol_copy.torsions):
    dihedral = simple_es_results.iloc[0][h+1]
    i,j,k,l = tor.indices
    RHS = tor.right_mask
    mol_copy.ase_molecule.set_dihedral(a1 = i,
                                       a2 = j,
                                       a3 = k,
                                       a4 = l,
                                       angle= float(dihedral),
                                       mask=RHS)
mol_copy.ase_molecule.set_calculator(Gaussian())

opt = optimize.BFGS(atoms=mol_copy.ase_molecule)
opt.run(fmax=0.05, trajectory='./drug_conformer_results/simple_es_paclitaxel.traj')
mol_copy.update_geometry_from_ase_mol()
simple_es_pickle = open("./drug_conformer_results/simple_es_final.pkl")
pickle.dump(mol_copy, simple_es_pickle)
