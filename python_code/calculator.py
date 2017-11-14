import rdkit, rdkit.Chem.rdDistGeom, rdkit.DistanceGeometry

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit import rdBase
from rdkit.Chem.rdMolTransforms import *


# Molecule calculators
def calc_energy(mol):
    """
    A function designed to calculate the energy of a molecule or TS.

    Input:
    * rdkit molecule object

    Output:
    * energy of the geometry calculated using rdkit's forcefield method
    * the resultant rdkit molecule object
    """
    # Setting the force field parameters
    try:
        molprop = rdkit.Chem.ChemicalForceFields.MMFFGetMoleculeProperties(mol)
    except:
        # Picking an arbitrarly large molecule to create a dummy mol properties forcefield object
        # This is primarily used for TSs
        test_mol = Chem.AddHs(Chem.MolFromSmiles(
            "COC(=O)Cc1ccccc1CC2CC(O)C(C(=C)CC#CCC=O)C2OO"
        ))
        molprop = rdkit.Chem.ChemicalForceFields.MMFFGetMoleculeProperties(test_mol)

    ff = rdkit.Chem.ChemicalForceFields.MMFFGetMoleculeForceField(mol, molprop)



    return ff.CalcEnergy(), mol

def constrained_energy_calc(mol, list_of_torsions, angles):
    """
    A function designed to calculate the energy of a molecule or TS.
    For this function, the molecule torsions are fixed,
    but the rest of the molecule is allowed to relax to a local minimum

    Input:
    * rdkit molecule object
    * list of torsions as (i,j,k,l) tuples
    * list of angles to set the corresponding torsions to

    Output:
    * energy of the geometry calculated using rdkit's forcefield method
    * the resultant rdkit molecule object
    """
    e, m = calc_energy(mol)
    # Setting the force field parameters
    try:
        molprop = rdkit.Chem.ChemicalForceFields.MMFFGetMoleculeProperties(mol)
    except:

        # Picking an arbitrarly large molecule to create a dummy mol properties forcefield object
        # This is primarily used for TSs
        test_mol = Chem.AddHs(Chem.MolFromSmiles(
            "COC(=O)Cc1ccccc1CC2CC(O)C(C(=C)CC#CCC=O)C2OO"
        ))
        molprop = rdkit.Chem.ChemicalForceFields.MMFFGetMoleculeProperties(test_mol)
    ff = rdkit.Chem.ChemicalForceFields.MMFFGetMoleculeForceField(mol, molprop)

    # Zipping the torsion list and the angles together
    geometries = zip(list_of_torsions, angles)

    # Generating a conformer to edit
    tup = mol.GetConformers()
    conformer = tup[0]

    # Setting the corresponding torsion to their proper angles
    for geometry in geometries:

        i = geometry[0][0]
        j = geometry[0][1]
        k = geometry[0][2]
        l = geometry[0][3]
        angle = geometry[1]

        SetDihedralDeg(conformer,
                           i,
                           j,
                           k,
                           l,
                           angle)

        # Constraining the torsions of the molecule
        ff.MMFFAddTorsionConstraint(i,j,k,l, True, angle, angle, 1)


    #e, m = calc_energy(mol)
    # Optimizing the molecule
    AllChem.MMFFOptimizeMolecule(mol)

    return ff.CalcEnergy(), mol

def optimize_and_calc(mol, list_of_torsions, angles):
    """
    A function designed to calculate the energy of a molecule or TS.
    For this function, the torsions of the molecule are set to their corresponding angles,
    then the molecule is allowed to relax to a local minimum

    Input:
    * rdkit molecule object
    * list of torsions as (i,j,k,l) tuples
    * list of angles to set the corresponding torsions to

    Output:
    * energy of the geometry calculated using rdkit's forcefield method
    * the resultant rdkit molecule object
    """

    # Setting the force field parameters
    try:
        molprop = rdkit.Chem.ChemicalForceFields.MMFFGetMoleculeProperties(mol)
    except:
        # Picking an arbitrarly large molecule to create a dummy mol properties forcefield object
        # This is primarily used for TSs
        test_mol = Chem.AddHs(Chem.MolFromSmiles(
            "COC(=O)Cc1ccccc1CC2CC(O)C(C(=C)CC#CCC=O)C2OO"
        ))

        molprop = rdkit.Chem.ChemicalForceFields.MMFFGetMoleculeProperties(test_mol)
    ff = rdkit.Chem.ChemicalForceFields.MMFFGetMoleculeForceField(mol, molprop)

    # Zipping the torsions and their corresponding angles together
    geometries = zip(list_of_torsions, angles)

    # Generating conformers and selecting one
    tup = mol.GetConformers()
    conformer = tup[0]

    for geometry in geometries:
        # Setting the torsion angles
        i = geometry[0][0]
        j = geometry[0][1]
        k = geometry[0][2]
        l = geometry[0][3]
        angle = geometry[1]

        SetDihedralDeg(conformer,
                           i,
                           j,
                           k,
                           l,
                           angle)

    # Allowing the molecule to relax
    AllChem.MMFFOptimizeMolecule(mol)

    return ff.CalcEnergy(), mol

# ~~~~~~~~~~~~~~~~~~~~~~

# TS calculators

def calc_energy_ts(mol):
    """
    A function designed to calculate the energy of a molecule or TS.

    Input:
    * rdkit molecule object

    Output:
    * energy of the geometry calculated using rdkit's forcefield method
    * the resultant rdkit molecule object
    """
    ff = rdkit.Chem.ChemicalForceFields.UFFGetMoleculeForceField(mol)

    return ff.CalcEnergy(), mol

def constrained_energy_calc_ts(mol, list_of_torsions, angles):
    """
    A function designed to calculate the energy of a molecule or TS.
    For this function, the molecule torsions are fixed,
    but the rest of the molecule is allowed to relax to a local minimum

    Input:
    * rdkit molecule object
    * list of torsions as (i,j,k,l) tuples
    * list of angles to set the corresponding torsions to

    Output:
    * energy of the geometry calculated using rdkit's forcefield method
    * the resultant rdkit molecule object
    """
    ff = rdkit.Chem.ChemicalForceFields.UFFGetMoleculeForceField(mol)

    # Zipping the torsion list and the angles together
    geometries = zip(list_of_torsions, angles)

    # Generating a conformer to edit
    tup = mol.GetConformers()
    conformer = tup[0]

    # Setting the corresponding torsion to their proper angles
    for geometry in geometries:

        i = geometry[0][0]
        j = geometry[0][1]
        k = geometry[0][2]
        l = geometry[0][3]
        angle = geometry[1]

        SetDihedralDeg(conformer,
                           i,
                           j,
                           k,
                           l,
                           angle)

        # Constraining the torsions of the molecule
        ff.UFFAddTorsionConstraint(i,j,k,l, True, angle, angle, 1)

    #e, m = calc_energy(mol)
    # Optimizing the molecule
    AllChem.UFFOptimizeMolecule(mol)

    return ff.CalcEnergy(), mol

def optimize_and_calc_ts(mol, list_of_torsions, angles):
    """
    A function designed to calculate the energy of a molecule or TS.
    For this function, the torsions of the molecule are set to their corresponding angles,
    then the molecule is allowed to relax to a local minimum

    Input:
    * rdkit molecule object
    * list of torsions as (i,j,k,l) tuples
    * list of angles to set the corresponding torsions to

    Output:
    * energy of the geometry calculated using rdkit's forcefield method
    * the resultant rdkit molecule object
    """

    ff = rdkit.Chem.ChemicalForceFields.UFFGetMoleculeForceField(mol)

    # Zipping the torsions and their corresponding angles together
    geometries = zip(list_of_torsions, angles)

    # Generating conformers and selecting one
    tup = mol.GetConformers()
    conformer = tup[0]

    for geometry in geometries:
        # Setting the torsion angles
        i = geometry[0][0]
        j = geometry[0][1]
        k = geometry[0][2]
        l = geometry[0][3]
        angle = geometry[1]

        SetDihedralDeg(conformer,
                           i,
                           j,
                           k,
                           l,
                           angle)

    # Allowing the molecule to relax
    AllChem.UFFOptimizeMolecule(mol)

    return ff.CalcEnergy(), mol
