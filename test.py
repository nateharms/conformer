def get_ts_rdmol(input_reaction):
    """

    """
    rxnFamilies = ["H_Abstraction"]
    rSpecies1, rSpecies2 = input_reaction.reactants
    pSpecies1, pSpecies2 = input_reaction.products

    for species in itertools.chain(input_reaction.reactants, input_reaction.products):
            species = species.generateResonanceIsomers()

    testReaction = Reaction(
        reactants = input_reaction.reactants,
        products = input_reaction.products,
        reversible = True)

    reactants = [species for species in input_reaction.reactants]
    products = [species for species in input_reaction.products]

    reactionList = []

    checkRxn = rmgDatabase.kinetics.generateReactionsFromFamilies(
        reactants,
        products,
        only_families=rxnFamilies)

    for rxn0 in checkRxn:
        reactionList.append(rxn0)
    logging.info("generateReactionsFromFamilies successfuly performed")

    gotOne=False
    for reaction in reactionList:
        # Check if any of the RMG proposed reactions matches the reaction in the mechanism
        if testReaction.isIsomorphic(reaction):
            print "Found matching reaction"
            # Now add the labeled atoms to the Molecule, and check all labels were added
            atLblsR = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])
            atLblsP = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])

            for reactant in reaction.reactants:
                reactant.clearLabeledAtoms()
                for atom in reactant.atoms:
                    for atomLabel in reaction.labeledAtoms:
                        if atom==atomLabel[1]:
                            atom.label = atomLabel[0]
                            atLblsR[atomLabel[0]] = True

            for product in reaction.products:
                product.clearLabeledAtoms()
                for atom in product.atoms:
                    for atomLabel in reaction.labeledAtoms:
                        if atom==atomLabel[1]:
                            atom.label = atomLabel[0]
                            atLblsP[atomLabel[0]] = True

            if all( atLblsR.values() ) and all( atLblsP.values() ):
                # We successfully labeled all of the atoms
                gotOne=True
                break

    rxn = QMReaction(reaction=reaction, settings=settings, tsDatabase=tsDatabase)


    reactant, product = rxn.setupMolecules()
    mol = reactant.toRDKitMol(removeHs=False)

    AllChem.EmbedMolecule(mol)
    labels, atomMatch = rxn.getLabels(reactant)
    tsRDMol, bm, rxn.reactantGeom = rxn.generateBoundsMatrix(reactant)
    bm = rxn.editMatrix(reactant, bm, labels)
    tsRDMol = rxn.reactantGeom.rd_embed(tsRDMol, 15, bm=bm, match=atomMatch)[0]
    tsRDMol = rdkit.Chem.rdchem.RWMol(tsRDMol)

    return tsRDMol

# ~~~~~~~~~~~~

def creat_pseudo_geometry(tsRDMol):
    tsRDMol_copy = tsRDMol.__copy__()

    for atom in tsRDMol_copy.GetAtoms():
        idx = atom.GetIdx()
        num = atom.GetAtomicNum()
        rmg_atom = reactant.atoms[idx]

        if rmg_atom.label:
            if rmg_atom.label == "*1":
                atom1_star = atom
            if rmg_atom.label == "*2":
                atom2_star = atom
            if rmg_atom.label == "*3":
                atom3_star = atom

    bond_between_23 = False
    try:
        tsRDMol_copy.AddBond(atom1_star.GetIdx(), atom2_star.GetIdx(), order=rdkit.Chem.rdchem.BondType.SINGLE)
    except RuntimeError:
        print "Bond already exists betwee 1* and 2*"
        bond_between_23 = True
        tsRDMol_copy.AddBond(atom2_star.GetIdx(), atom3_star.GetIdx(), order=rdkit.Chem.rdchem.BondType.SINGLE)

    return tsRDMol_copy

#~~~~~~~~~


def get_ts_torsion_energies(tsRDMol_copy, delta):

get_torsion_list(tsRDMol_copy)





torsion_combos = list( itertools.combinations_with_replacement(torsion_angles, len(torsion_list)) )
    if len(torsion_list) != 1:
        torsion_combos = list(
            set(
                torsion_combos +
                list(itertools.combinations_with_replacement(
                    torsion_angles[::-1], len(torsion_list)
                ))))

tup = tsRDMol_copy.GetConformers()
conformer = tup[0]

# Picking a random torsion combo
for combo in torsion_combos:
    geometry = zip(torsion_list, combo)
    for torsion in geometry:
        print torsion
        i = torsion[0][0]
        j = torsion[0][1]
        k = torsion[0][2]
        l = torsion[0][3]
        angle = torsion[1]

        SetDihedralDeg(conformer,
                       i,
                       j,
                       k,
                       l,
                       angle)

    if bond_between_23 == False:
        tsRDMol_copy.RemoveBond(atom1_star.GetIdx(), atom2_star.GetIdx())
    else:
        tsRDMol_copy.RemoveBond(atom2_star.GetIdx(), atom3_star.GetIdx())

    print calc_energy(tsRDMol_copy)
