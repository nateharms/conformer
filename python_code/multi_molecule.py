import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolTransforms import *
import ase
from ase import Atom, Atoms
import rmgpy
from rmgpy.molecule import Molecule
import py3Dmol
import numpy as np
from numpy import pi

class Torsion():
    """
    A class that acts as a container for torsion information.
    The information stored is as follows:

    * indicies (list of ints): indicies of atoms that describe the dihedral
    * dihedral (float): the angle (in degrees) of the dihedral
    * LHS (list of ints): indices of atoms that are branced off the left side of the torsion
    * RHS (list of ints): indices of atoms that are branced off the right side of the torsion
    """

    def __init__(self, indices, dihedral, LHS, RHS):
        self.indices = indices
        self.dihedral = dihedral
        self.LHS = LHS
        self.RHS = RHS




class Multi_Molecule():
    """
    A class that allows for one to create RMG, RDKit and ASE
    molecules from a single string with identical atom indicies

    Inputs:
    * smiles (str): a SMILES string that describes the molecule of interest
    """
    def __init__(self, smiles):

        self.smiles = smiles

        # Make these more explicit (make rmg mol from smiles)

        self.get_rmg_molecule()
        self.get_rdkit_molecule()
        self.set_rmg_coords("RDKit")
        self.get_ase_molecule()
        self.get_torsion_list()
        self.get_torsions()


    def get_rmg_molecule(self):
        """
        A method to obtain and set the RMG Molecule
        """
        self.rmg_molecule = Molecule(SMILES=self.smiles)
    def get_rdkit_molecule(self):
        """
        A method to create an RDKit Molecule from the rmg_molecule.
        Indicies will be the same as in the RMG Molecule
        """

        RDMol = self.rmg_molecule.toRDKitMol(removeHs=False)

        rdkit.Chem.AllChem.EmbedMolecule(RDMol)

        self.rdkit_molecule = RDMol

    def get_ase_molecule(self):
        """
        A method to create an ASE Molecule from the rdkit_molecule.
        Indicies will be the same as in the RMG and RDKit Molecule.
        """
        mol_list = AllChem.MolToMolBlock(self.rdkit_molecule).split('\n')
        ase_atoms = []
        for i, line in enumerate(mol_list):

            if i > 3:

                try:
                    atom0, atom1, bond, rest = line
                    atom0 = int(atom0)
                    atom0 = int(atom1)
                    bond = float(bond)

                except ValueError:
                    try:
                        x, y, z, symbol = line.split()[0:4]
                        x = float(x)
                        y = float(y)
                        z = float(z)
                        #print symbol

                        ase_atoms.append(Atom(symbol=symbol, position=(x,y,z)))

                    except:
                        continue

        self.ase_molecule = Atoms(ase_atoms)
        return self.ase_molecule

    def view_mol(self):
        """
        A method designed to create a 3D figure of the Multi_Molecule with py3Dmol
        """
        mb  = Chem.MolToMolBlock(self.rdkit_molecule)
        p = py3Dmol.view(width=400, height=400)
        p.addModel(mb, "sdf")
        p.setStyle({'stick':{}})
        p.setBackgroundColor('0xeeeeee')
        p.zoomTo()
        return p.show()


    def get_torsion_list(self):
        """
        A method to return a list of the possible torsions in a Multi_Molecule.
        This uses the RDKit framework to do this.
        """
        RDMol = self.rdkit_molecule
        torsion_list = []
        for bond1 in RDMol.GetBonds():
            atom1 = bond1.GetBeginAtom()
            atom2 = bond1.GetEndAtom()
            if atom1.IsInRing() or atom2.IsInRing():
                # Making sure that bond1 we're looking at are in a ring
                continue

            bond_list1 = list(atom1.GetBonds())
            bond_list2 = list(atom2.GetBonds())

            if not len(bond_list1) > 1 and not len(bond_list2) > 1:
                # Making sure that there are more than one bond attached to
                # the atoms we're looking at
                continue

            # Getting the 0th and 3rd atom and insuring that atoms
            # attached to the 1st and 2nd atom are not terminal hydrogens
            # We also make sure that all of the atoms are properly bound together

            # If the above are satisified, we append a tuple of the torsion our torsion_list
            got_atom0 = False
            got_atom3 = False

            for bond0 in bond_list1:
                atomX = bond0.GetOtherAtom(atom1)
                if atomX.GetAtomicNum() == 1 and len(atomX.GetBonds()) == 1:
                    # This means that we have a terminal hydrogen, skip this
                    # NOTE: for H_abstraction TSs, a non teminal H should exist
                    continue
                if atomX.GetIdx() != atom2.GetIdx():
                    got_atom0 = True
                    atom0 = atomX

            for bond2 in bond_list2:
                atomY = bond2.GetOtherAtom(atom2)
                if atomY.GetAtomicNum() == 1 and len(atomY.GetBonds()) == 1:
                    # This means that we have a terminal hydrogen, skip this
                    continue
                if atomY.GetIdx() != atom1.GetIdx():
                    got_atom3 = True
                    atom3 = atomY

            if not (got_atom0 and got_atom3):
                # Making sure atom0 and atom3 were not found
                continue

            # Looking to make sure that all of the atoms are properly bonded to eached
            if (
                RDMol.GetBondBetweenAtoms(atom0.GetIdx(), atom1.GetIdx()) and
                RDMol.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx()) and
                RDMol.GetBondBetweenAtoms(atom2.GetIdx(), atom3.GetIdx())   ) :

                torsion_tup = (atom0.GetIdx(), atom1.GetIdx(), atom2.GetIdx(), atom3.GetIdx())
                torsion_list.append(torsion_tup)

        self.torsion_list = torsion_list
        return self.torsion_list

    def get_torsions(self):
        torsions = []
        for indices in self.torsion_list:
            i, j, k, l = indices

            dihedral = self.ase_molecule.get_dihedral(i,j,k,l)
            tor = Torsion(indices=indices, dihedral=dihedral, LHS=[], RHS=[])
            LHS = self.get_LHS(tor)
            RHS = self.get_RHS(tor)

            torsions.append(Torsion(indices, dihedral, LHS, RHS))
        self.torsions = torsions
        return self.torsions

    def get_LHS(self, Torsion):

        rdkit_atoms = self.rdkit_molecule.GetAtoms()

        L1, L0, R0, R1 = Torsion.indices
        rd_atom_L1 = rdkit_atoms[L1]
        rd_atom_L0 = rdkit_atoms[L0]
        rd_atom_R0 = rdkit_atoms[R0]
        rd_atom_R1 = rdkit_atoms[R1]

        # trying to get the left hand side of this torsion
        LHS_atoms_index = [rd_atom_L0.GetIdx(), rd_atom_L1.GetIdx()]
        RHS_atoms_index = [rd_atom_R0.GetIdx(), rd_atom_R1.GetIdx()]

        complete_LHS = False
        i = 0
        atom_index = LHS_atoms_index[0]
        while complete_LHS == False:
            try:
                LHS_atom = rdkit_atoms[atom_index]
                for neighbor in LHS_atom.GetNeighbors():
                    if (neighbor.GetIdx() in LHS_atoms_index) or (neighbor.GetIdx() in RHS_atoms_index):
                        continue
                    else:
                        LHS_atoms_index.append(neighbor.GetIdx())
                i +=1
                atom_index = LHS_atoms_index[i]

            except IndexError:
                complete_LHS = True

        return LHS_atoms_index

    def get_RHS(self, Torsion):

        rdkit_atoms = self.rdkit_molecule.GetAtoms()

        L1, L0, R0, R1 = Torsion.indices
        rd_atom_L1 = rdkit_atoms[L1]
        rd_atom_L0 = rdkit_atoms[L0]
        rd_atom_R0 = rdkit_atoms[R0]
        rd_atom_R1 = rdkit_atoms[R1]

        # trying to get the left hand side of this torsion
        LHS_atoms_index = [rd_atom_L0.GetIdx(), rd_atom_L1.GetIdx()]
        RHS_atoms_index = [rd_atom_R0.GetIdx(), rd_atom_R1.GetIdx()]

        complete_RHS = False
        i = 0
        atom_index = RHS_atoms_index[0]
        while complete_RHS == False:
            try:
                RHS_atom = rdkit_atoms[atom_index]
                for neighbor in RHS_atom.GetNeighbors():
                    if (neighbor.GetIdx() in RHS_atoms_index) or (neighbor.GetIdx() in LHS_atoms_index):
                        continue
                    else:
                        RHS_atoms_index.append(neighbor.GetIdx())
                i +=1
                atom_index = RHS_atoms_index[i]

            except IndexError:
                complete_RHS = True

        return RHS_atoms_index

    def set_rmg_coords(self, molecule_base):

        if molecule_base == "RDKit":
            mol_list = AllChem.MolToMolBlock(self.rdkit_molecule).split('\n')
            for i, atom in enumerate(self.rmg_molecule.atoms):

                j = i + 4
                coords = mol_list[j].split()[:3]

                for k, coord in enumerate(coords):

                    coords[k] = float(coord)
                atom.coords = np.array(coords)

        elif molecule_base == "ASE":
            for i, position in enumerate(self.ase_molecule.get_positions()):
                self.rmg_molecule.atoms[i].coords = position




    def update_geometry_from_rdkit_mol(self):
        # In order to update the ase molecule you simply need to rerun the get_ase_molecule method
        self.get_ase_molecule()
        self.set_rmg_coords("RDKit")



    def update_geometry_from_ase_mol(self):

        self.set_rmg_coords("ASE")

        #setting the geometries of the rdkit molecule

        positions = self.ase_molecule.get_positions()

        conf = self.rdkit_molecule.GetConformers()[0]

        for i, atom in enumerate(self.rdkit_molecule.GetAtoms()):
            conf.SetAtomPosition(i, positions[i])


    def update_geometry_from_rmg_mol(self):
        # I don't know why you would ever want to do this, but okay...

        conf = self.rdkit_molecule.GetConformers()[0]
        ase_atoms = []
        for i, atom in enumerate(self.rmg_molecule.atoms):
            x,y,z = atom.coords
            symbol = atom.symbol

            conf.SetAtomPosition(i, [x,y,z])


            ase_atoms.append(Atom(symbol=symbol, position=(x, y, z)))

        self.ase_molecule = Atoms(ase_atoms)
