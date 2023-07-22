"""
Module for generating rdkit molobj/smiles/molecular graph from free atoms

Implementation by Jan H. Jensen, based on the paper

    Yeonjoon Kim and Woo Youn Kim
    "Universal Structure Conversion Method for Organic Molecules: From Atomic Connectivity
    to Three-Dimensional Geometry"
    Bull. Korean Chem. Soc. 2015, Vol. 36, 1769-1777
    DOI: 10.1002/bkcs.10334

"""

import copy
import itertools

from rdkit.Chem import rdmolops
from rdkit.Chem import rdchem
try:
    from rdkit.Chem import rdEHTTools #requires RDKit 2019.9.1 or later
except ImportError:
    rdEHTTools = None
    
from collections import defaultdict, Counter
from itertools import combinations

import numpy as np
import networkx as nx

from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
import sys

global __ATOM_LIST__
__ATOM_LIST__ = \
    ['h',  'he',
     'li', 'be', 'b',  'c',  'n',  'o',  'f',  'ne',
     'na', 'mg', 'al', 'si', 'p',  's',  'cl', 'ar',
     'k',  'ca', 'sc', 'ti', 'v ', 'cr', 'mn', 'fe', 'co', 'ni', 'cu',
     'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr',
     'rb', 'sr', 'y',  'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag',
     'cd', 'in', 'sn', 'sb', 'te', 'i',  'xe',
     'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy',
     'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w',  're', 'os', 'ir', 'pt',
     'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn',
     'fr', 'ra', 'ac', 'th', 'pa', 'u',  'np', 'pu']


global atomic_valence
global atomic_valence_electrons
global atomic_electronegativity

atomic_valence = defaultdict(list)
atomic_valence[1] = [1]
atomic_valence[5] = [3,4]
atomic_valence[6] = [4]
atomic_valence[7] = [3,4]
atomic_valence[8] = [2]
atomic_valence[9] = [1]
atomic_valence[14] = [4]
atomic_valence[15] = [5,3] #[5,4,3]
atomic_valence[16] = [6,3,2] #[6,4,2]
atomic_valence[17] = [1]
atomic_valence[32] = [4]
atomic_valence[35] = [1]
atomic_valence[53] = [1]

atomic_valence_electrons = {}
atomic_valence_electrons[1] = 1
atomic_valence_electrons[5] = 3
atomic_valence_electrons[6] = 4
atomic_valence_electrons[7] = 5
atomic_valence_electrons[8] = 6
atomic_valence_electrons[9] = 7
atomic_valence_electrons[14] = 4
atomic_valence_electrons[15] = 5
atomic_valence_electrons[16] = 6
atomic_valence_electrons[17] = 7
atomic_valence_electrons[32] = 4
atomic_valence_electrons[35] = 7
atomic_valence_electrons[53] = 7

max_allowed_bonds = defaultdict(list)
max_allowed_bonds[1] = [1]
max_allowed_bonds[5] = [4]
max_allowed_bonds[6] = [4]
max_allowed_bonds[7] = [4]
max_allowed_bonds[8] = [3]
max_allowed_bonds[9] = [1]
max_allowed_bonds[14] = [4]
max_allowed_bonds[15] = [5] #[5,4,3]
max_allowed_bonds[16] = [6] #[6,4,2]
max_allowed_bonds[17] = [1]
max_allowed_bonds[32] = [4]
max_allowed_bonds[35] = [1]
max_allowed_bonds[53] = [1]

atomic_electronegativity = {}
atomic_electronegativity[1] = 2.2
atomic_electronegativity[5] = 2.04
atomic_electronegativity[6] = 2.55
atomic_electronegativity[7] = 3.04
atomic_electronegativity[8] = 3.44
atomic_electronegativity[9] = 3.98
atomic_electronegativity[14] = 1.90
atomic_electronegativity[15] = 2.19
atomic_electronegativity[16] = 2.58
atomic_electronegativity[17] = 3.16
atomic_electronegativity[35] = 2.96
atomic_electronegativity[53] = 2.66


def mol_with_atom_index(mol):
    atoms = mol.GetNumAtoms()
    for idx in range(atoms):
        mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx()))
    return mol

def str_atom(atom):
    """
    convert integer atom to string atom
    """
    global __ATOM_LIST__
    atom = __ATOM_LIST__[atom - 1]
    return atom


def int_atom(atom):
    """
    convert str atom to integer atom
    """
    global __ATOM_LIST__
    #print(atom)
    atom = atom.lower()
    return __ATOM_LIST__.index(atom) + 1


def get_UA(maxValence_list, valence_list, atomic_valence_electrons, atom_electroNegativiy): #we have to return the version where O has to lowest number of bonds (rather 2 than 3)
    """
    """
    UA = []
    DU = []
    octet_offset = []
    electro_negativity = []
    for i, (maxValence, valence, atom_val, EN) in enumerate(zip(maxValence_list, valence_list, atomic_valence_electrons, atom_electroNegativiy)):
        if not maxValence - valence > 0:
            continue
        UA.append(i)
        DU.append(maxValence - valence)
        octet_offset.append(atom_val-8+valence)
        electro_negativity.append(EN)
    return UA, DU, octet_offset, electro_negativity



def find_combinations(UA, UA_pairs, DU, octet_offset, electro_negativity): # need to adapt the feeding function to also submit the cases with only one ulti-bond to determine if double or tripl
    results = []
    DU_dict = dict(zip(UA, DU))
    octet_offset_dict = dict(zip(UA, octet_offset))
    electro_negativity_dict = dict(zip(UA, electro_negativity))

    # Create a set of keys that appear in UA_pair
    keys_in_UA_pair = set()
    for pair in UA_pairs:
        keys_in_UA_pair.update(pair)

    # Remove keys not in UA_pair from DU_dict
    for key in list(DU_dict.keys()):
        if key not in keys_in_UA_pair:
            del DU_dict[key]

    # Remove keys not in UA_pair from octet_offset_dict
    for key in list(octet_offset_dict.keys()):
        if key not in keys_in_UA_pair:
            del octet_offset_dict[key]

    # Remove keys not in UA_pair from electro_negativity_dict
    for key in list(electro_negativity_dict.keys()):
        if key not in keys_in_UA_pair:
            del electro_negativity_dict[key]

    #new collections
    overallMin_octet_nonZero_offset_atoms_count = float('inf')
    overallMin_octet_neg_offset_atoms_count = float('inf')
    overallMin_octet_pos_offset_atoms_count = float('inf')

    for r in range(1, len(UA_pairs) + 1):
        for UA_pair_subset in combinations(UA_pairs, r):
            flat_subset = [i for sub in UA_pair_subset for i in sub]
            counter = Counter(flat_subset)
            UA_pair_subset = {key: 1 for key in UA_pair_subset}

            #check if atoms are still unsaturated or oversaturated which in this case means how many more bonds we would have
            # compared to the maximum number of allowed bonds (e.g. 4 for N) (positive values mean unsaturated)
            new_UA_or_OA = {}
            for key, value in counter.items():
                diff = DU_dict.get(key, 0) - value
                new_UA_or_OA[key] = diff

            #at this point we dont allow more bonds the maximum allowed number set in beginning (in dict on top of this script)
            any_oversaturated = any(value < 0 for value in new_UA_or_OA.values())
            if any_oversaturated:
                continue

            #check how much we are off from the octet rule by atom (positive values means we have more bonds than we would need to satisfy octet rule)
            new_octet_offset = {}
            for key in octet_offset_dict:
                value_octet_offset = octet_offset_dict[key]
                value_counter = counter.get(key, 0)
                new_value = value_octet_offset + value_counter
                new_octet_offset[key] = new_value

            #determine potential triple bonds
            for UA_tuple in UA_pair_subset.keys():
                keys_in_tuple = UA_tuple
                negative_values = [new_octet_offset.get(key, 0) for key in keys_in_tuple if
                                   new_octet_offset.get(key, 0) < 0]
                if len(negative_values) == 2:
                    UA_pair_subset[UA_tuple] = 2 #triple bonds will be -1

            for key in new_octet_offset:
                for ua_key in UA_pair_subset:
                    if key in ua_key and UA_pair_subset[ua_key] == 2:
                        new_octet_offset[key] += UA_pair_subset[ua_key]

            #check if the atoms with the highest electronegativity have the electrons
            summed_electroneg_UA = 0 #summed electronegativity for unsaturated atoms
            summed_electroneg_SA = 0 #summed electronegativity for saturated atoms

            for key in new_octet_offset:
                if new_octet_offset[key] < 0:
                    summed_electroneg_UA += electro_negativity_dict[key]
                else:
                    summed_electroneg_SA += electro_negativity_dict[key]

            summed_electroneg_UA_SA_diff = summed_electroneg_UA - summed_electroneg_SA #when negative the more eletronegative atoms have the electrons




            #check for number of atoms that dont satisfy octet and the highest offset. positive values mean that we have more bonds than needed for
            octet_nonZero_offset_atoms_count = 0
            octet_neg_offset_atoms_count = 0 #missing bonds
            highest_neg_octet_offset = float('inf') #will result in positive or negative charges
            octet_pos_offset_atoms_count = 0  #too many bonds
            highest_pos_octet_offset = float('-inf') #will result in positive charges
            for value in new_octet_offset.values():
                if value != 0:
                    octet_nonZero_offset_atoms_count += 1
                if value < 0:
                    octet_neg_offset_atoms_count += 1
                if value < highest_neg_octet_offset:
                    highest_neg_octet_offset = value
                if value > 0:
                    octet_pos_offset_atoms_count += 1
                if value > highest_pos_octet_offset:
                    highest_pos_octet_offset = value


            overallMin_octet_nonZero_offset_atoms_count = min(overallMin_octet_nonZero_offset_atoms_count, octet_nonZero_offset_atoms_count)


            # Store the difference together with the subset
            results.append((octet_nonZero_offset_atoms_count, summed_electroneg_UA_SA_diff, UA_pair_subset))

    #make sure best octet rule saturation is used
    UA_pairs_satisfying_maxBonds_and_minOctetViolations = [[t] for t in results if t[0] == overallMin_octet_nonZero_offset_atoms_count]


    #select the combination where the electronegtivity determines where the electrons should be
    best_electroneg_diff = min(range(len(UA_pairs_satisfying_maxBonds_and_minOctetViolations)),
                             key=lambda index: UA_pairs_satisfying_maxBonds_and_minOctetViolations[index][0][1])

    UA_pairs_satisfying_maxBonds_and_minOctetViolations = [UA_pairs_satisfying_maxBonds_and_minOctetViolations[best_electroneg_diff][0][2]]

    return UA_pairs_satisfying_maxBonds_and_minOctetViolations

def get_BO(AC, UA, DU, valences, UA_pairs, use_graph=True):
    """
    """
    BO = AC.copy()

    for (i, j), value in UA_pairs.items():
        BO[i, j] += value
        BO[j, i] += value
        pass

    return BO


def valences_not_too_large(BO, valences):
    """
    """
    number_of_bonds_list = BO.sum(axis=1)
    for valence, number_of_bonds in zip(valences, number_of_bonds_list):
        if number_of_bonds > valence:
            return False

    return True

def charge_is_OK(BO, AC, charge, DU, atomic_valence_electrons, atoms, valences,
                 allow_charged_fragments=True):
    # total charge
    Q = 0

    # charge fragment list
    q_list = []

    if allow_charged_fragments:

        BO_valences = list(BO.sum(axis=1))
        for i, atom in enumerate(atoms):
            q, radical = get_atomic_charge(atom, atomic_valence_electrons[atom], BO_valences[i])
            Q += q
            if atom == 6:
                number_of_single_bonds_to_C = list(BO[i, :]).count(1)
                if number_of_single_bonds_to_C == 2 and BO_valences[i] == 2:
                    Q += 1
                    q = 2
                if number_of_single_bonds_to_C == 3 and Q + 1 < charge:
                    Q += 2
                    q = 1

            if q != 0:
                q_list.append(q)

    return (charge == Q)

def BO_is_OK(BO, AC, charge, DU, atomic_valence_electrons, atoms, valences,
    allow_charged_fragments=True):
    """
    Sanity of bond-orders

    args:
        BO -
        AC -
        charge -
        DU - 


    optional
        allow_charges_fragments - 


    returns:
        boolean - true of molecule is OK, false if not
    """

    if not valences_not_too_large(BO, valences):
        return False

    check_sum = (BO - AC).sum() == sum(DU)
    check_charge = charge_is_OK(BO, AC, charge, DU, atomic_valence_electrons, atoms, valences,
                                allow_charged_fragments)

    if check_charge and check_sum: 
        return True

    return False


def get_atomic_charge(atom, atomic_valence_electrons, BO_valence):
    """
    """
    radical = False

    if atom == 1:
        charge = 1 - BO_valence
    elif atom == 5:
        charge = 3 - BO_valence
    elif atom == 6 and BO_valence == 2:
        charge = 0
    elif atom == 7 and BO_valence == 2:
        charge = 1
        radical = True
    elif atom == 15 and BO_valence == 5:
        charge = 0
    elif atom == 16 and BO_valence == 6:
        charge = 0
    else:
        charge = atomic_valence_electrons - 8 + BO_valence

    return charge, radical


def clean_charges(mol):
    """
    This hack should not be needed anymore, but is kept just in case

    """

    Chem.SanitizeMol(mol)
    #rxn_smarts = ['[N+:1]=[*:2]-[C-:3]>>[N+0:1]-[*:2]=[C-0:3]',
    #              '[N+:1]=[*:2]-[O-:3]>>[N+0:1]-[*:2]=[O-0:3]',
    #              '[N+:1]=[*:2]-[*:3]=[*:4]-[O-:5]>>[N+0:1]-[*:2]=[*:3]-[*:4]=[O-0:5]',
    #              '[#8:1]=[#6:2]([!-:6])[*:3]=[*:4][#6-:5]>>[*-:1][*:2]([*:6])=[*:3][*:4]=[*+0:5]',
    #              '[O:1]=[c:2][c-:3]>>[*-:1][*:2][*+0:3]',
    #              '[O:1]=[C:2][C-:3]>>[*-:1][*:2]=[*+0:3]']

    rxn_smarts = ['[#6,#7:1]1=[#6,#7:2][#6,#7:3]=[#6,#7:4][CX3-,NX3-:5][#6,#7:6]1=[#6,#7:7]>>'
                  '[#6,#7:1]1=[#6,#7:2][#6,#7:3]=[#6,#7:4][-0,-0:5]=[#6,#7:6]1[#6-,#7-:7]',
                  '[#6,#7:1]1=[#6,#7:2][#6,#7:3](=[#6,#7:4])[#6,#7:5]=[#6,#7:6][CX3-,NX3-:7]1>>'
                  '[#6,#7:1]1=[#6,#7:2][#6,#7:3]([#6-,#7-:4])=[#6,#7:5][#6,#7:6]=[-0,-0:7]1']

    fragments = Chem.GetMolFrags(mol,asMols=True,sanitizeFrags=False)

    for i, fragment in enumerate(fragments):
        for smarts in rxn_smarts:
            patt = Chem.MolFromSmarts(smarts.split(">>")[0])
            while fragment.HasSubstructMatch(patt):
                rxn = AllChem.ReactionFromSmarts(smarts)
                ps = rxn.RunReactants((fragment,))
                fragment = ps[0][0]
                Chem.SanitizeMol(fragment)
        if i == 0:
            mol = fragment
        else:
            mol = Chem.CombineMols(mol, fragment)

    return mol


def BO2mol(mol, BO_matrix, atoms, atomic_valence_electrons,
           mol_charge, allow_charged_fragments=True,  use_atom_maps=False):
    """
    based on code written by Paolo Toscani

    From bond order, atoms, valence structure and total charge, generate an
    rdkit molecule.

    args:
        mol - rdkit molecule
        BO_matrix - bond order matrix of molecule
        atoms - list of integer atomic symbols
        atomic_valence_electrons -
        mol_charge - total charge of molecule

    optional:
        allow_charged_fragments - bool - allow charged fragments

    returns
        mol - updated rdkit molecule with bond connectivity

    """

    l = len(BO_matrix)
    l2 = len(atoms)
    BO_valences = list(BO_matrix.sum(axis=1))

    if (l != l2):
        raise RuntimeError('sizes of adjMat ({0:d}) and Atoms {1:d} differ'.format(l, l2))

    rwMol = Chem.RWMol(mol)


    bondTypeDict = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE
    }

    for i in range(l):
        for j in range(i + 1, l):
            bo = int(round(BO_matrix[i, j]))
            if (bo == 0):
                continue
            bt = bondTypeDict.get(bo, Chem.BondType.SINGLE)
            rwMol.AddBond(i, j, bt)

    mol = rwMol.GetMol()


    if allow_charged_fragments:
        mol = set_atomic_charges(
            mol,
            atoms,
            atomic_valence_electrons,
            BO_valences,
            BO_matrix,
            mol_charge,
            use_atom_maps)
    else:
        mol = set_atomic_radicals(mol, atoms, atomic_valence_electrons, BO_valences,
                                                            use_atom_maps)

    return mol


def set_atomic_charges(mol, atoms, atomic_valence_electrons,
                       BO_valences, BO_matrix, mol_charge,
                       use_atom_maps):
    """
    """
    q = 0
    for i, atom in enumerate(atoms):
        a = mol.GetAtomWithIdx(i)
        if use_atom_maps:
            a.SetAtomMapNum(i+1)
        charge, radical = get_atomic_charge(atom, atomic_valence_electrons[atom], BO_valences[i])

        q += charge
        if atom == 6:
            number_of_single_bonds_to_C = list(BO_matrix[i, :]).count(1)
            if number_of_single_bonds_to_C == 2 and BO_valences[i] == 2:
                q += 1
                charge = 0
            if number_of_single_bonds_to_C == 3 and q + 1 < mol_charge:
                q += 2
                charge = 1
            if BO_matrix[i, :].sum() == 3:
                charge = 1
            #if charge == -1:
            #    pass

        if (abs(charge) != 0 and radical == False):
            a.SetFormalCharge(int(charge))
        if (abs(charge) != 0 and radical == True):
            a.SetFormalCharge(int(charge))
            a.SetNumRadicalElectrons(int(atomic_valence_electrons[atom] - BO_valences[i] - charge))
    #mol = clean_charges(mol)

    return mol


def set_atomic_radicals(mol, atoms, atomic_valence_electrons, BO_valences,
                                                use_atom_maps):
    """

    The number of radical electrons = absolute atomic charge

    """
    for i, atom in enumerate(atoms):
        a = mol.GetAtomWithIdx(i)
        if use_atom_maps:
            a.SetAtomMapNum(i+1)
        charge = get_atomic_charge(
            atom,
            atomic_valence_electrons[atom],
            BO_valences[i])

        if (abs(charge) > 0):
            a.SetNumRadicalElectrons(abs(int(charge)))

    return mol


def get_bonds(UA, AC):
    """

    """
    bonds = []

    for k, i in enumerate(UA):
        for j in UA[k + 1:]:
            if AC[i, j] == 1:
                bonds.append(tuple(sorted([i, j])))

    return bonds


def get_UA_pairs(UA, AC, use_graph=True):
    """

    """

    bonds = get_bonds(UA, AC)

    if len(bonds) == 0:
        return [()]

    if use_graph:
        G = nx.Graph()
        G.add_edges_from(bonds)
        UA_pairs = [list(nx.max_weight_matching(G))]

        #Get conjugated atoms
        components = nx.connected_components(G)

        # Initialize an empty list to store the connections for each conjugated system
        component_edges = []

        # Iterate over the bonds
        for component in components:
            # Create a subgraph for the current component
            subgraph = G.subgraph(component)
            # Get the bonds of the subgraph and append to the list
            component_edges.append(list(subgraph.edges()))



        return component_edges

    max_atoms_in_combo = 0
    UA_pairs = [()]
    for combo in list(itertools.combinations(bonds, int(len(UA) / 2))):
        flat_list = [item for sublist in combo for item in sublist]
        atoms_in_combo = len(set(flat_list))
        if atoms_in_combo > max_atoms_in_combo:
            max_atoms_in_combo = atoms_in_combo
            UA_pairs = [combo]

        elif atoms_in_combo == max_atoms_in_combo:
            UA_pairs.append(combo)

    return UA_pairs

def is_nested_list(lst):
    return isinstance(lst, list) and any(isinstance(i, list) for i in lst)

def process_pairs(UA_pairs_list_raw, UA, DU_from_AC, octet_offset, electro_negativity):
    UA_pairs_list = []
    #for conjugated_bond_cluster in UA_pairs_list_raw:
    #    UA_pairs_list.extend(find_combinations(UA, conjugated_bond_cluster, DU_from_AC, octet_offset, electro_negativity))
    conjugated_bond_cluster = [item for sublist in UA_pairs_list_raw for item in sublist]
    pass
    if len(conjugated_bond_cluster) > 0:
        UA_pairs_list = find_combinations(UA, conjugated_bond_cluster, DU_from_AC, octet_offset, electro_negativity)
    else:
        return {}
    pass
    result_dict = {}
    for ua_dict in UA_pairs_list:
        result_dict.update(ua_dict)

    return result_dict

def AC2BO(AC, atoms, charge, allow_charged_fragments=True, use_graph=True):
    """

    implemenation of algorithm shown in Figure 2

    UA: unsaturated atoms

    DU: degree of unsaturation (u matrix in Figure)

    best_BO: Bcurr in Figure

    """
    AC_save = AC.copy()

    global atomic_valence
    global atomic_valence_electrons
    reps = 1

    attempt_to_keep_bonds = np.any(AC_save == 2)
    if attempt_to_keep_bonds == True:
        reps = 2

    remove_attempted_bonds = False
    valence_valid = False

    for try_con in range(reps):
        AC = AC_save.copy()
        if attempt_to_keep_bonds:
            if remove_attempted_bonds == False:
                AC[AC == 2] = 1
            if remove_attempted_bonds == True:
                AC[AC == 2] = 0

        # make a list of valences, e.g. for CO: [[4],[2,1]]
        valences_list_of_lists = []
        atom_valenceElectrons = []
        atom_electroNegativiy = []
        AC_valence = list(AC.sum(axis=1))

        for i,(atomicNum,valence) in enumerate(zip(atoms,AC_valence)):
            # valence can't be smaller than number of neighbourgs
            possible_valence = [x for x in max_allowed_bonds[atomicNum] if x >= valence]
            if not possible_valence:
                if attempt_to_keep_bonds == False:
                    print('Valence of atom',i,'is',valence,'which bigger than allowed max',max(atomic_valence[atomicNum]),'. Stopping')
                    #sys.exit()
                    return None
                if attempt_to_keep_bonds == True:
                    remove_attempted_bonds = True
                    break
            valences_list_of_lists.append(possible_valence)
            atom_valenceElectrons.append(atomic_valence_electrons[atomicNum])
            atom_electroNegativiy.append(atomic_electronegativity[atomicNum])

            if i == len(atoms)-1:
                valence_valid = True

        #trying replacement
        # convert [[4],[2,1]] to [[4,2],[4,1]]
        valences_list = itertools.product(*valences_list_of_lists)


        if valence_valid == True:
            break

    for valences in valences_list: #atm just 1 iteration because I use maximum allowed bonds
        UA, DU_from_AC, octet_offset, electro_negativity = get_UA(valences, AC_valence, atom_valenceElectrons, atom_electroNegativiy)

    UA_pairs_list_raw = get_UA_pairs(UA, AC, use_graph=use_graph)
    UA_pairs = process_pairs(UA_pairs_list_raw, UA, DU_from_AC, octet_offset, electro_negativity)

    BO = get_BO(AC, UA, DU_from_AC, valences, UA_pairs, use_graph=use_graph)


    return BO, atomic_valence_electrons, AC


def AC2mol(mol, AC, atoms, charge, allow_charged_fragments=True, 
           use_graph=True, use_atom_maps=False):
    """
    """

    # convert AC matrix to bond order (BO) matrix
    BO, atomic_valence_electrons, AC = AC2BO(
        AC,
        atoms,
        charge,
        allow_charged_fragments=allow_charged_fragments,
        use_graph=use_graph)

    # add BO connectivity and charge info to mol object
    mol = BO2mol(
        mol,
        BO,
        atoms,
        atomic_valence_electrons,
        charge,
        allow_charged_fragments=allow_charged_fragments,
        use_atom_maps=use_atom_maps)


    # BO2mol returns an arbitrary resonance form. Let's make the rest
    flags = rdchem.ResonanceFlags.ALLOW_INCOMPLETE_OCTETS | rdchem.ResonanceFlags.KEKULE_ALL
    mols = rdchem.ResonanceMolSupplier(mol, flags)

    #mols = rdchem.ResonanceMolSupplier(mol, Chem.UNCONSTRAINED_CATIONS, Chem.UNCONSTRAINED_ANIONS, ALLOW_INCOMPLETE_OCTETS = 1)

    mols = [mol for mol in mols]
    return mols, AC


def get_proto_mol(atoms):
    """
    """
    mol = Chem.MolFromSmarts("[#" + str(atoms[0]) + "]")
    rwMol = Chem.RWMol(mol)
    for i in range(1, len(atoms)):
        a = Chem.Atom(atoms[i])
        rwMol.AddAtom(a)

    mol = rwMol.GetMol()

    return mol


def read_xyz_file(filename, look_for_charge=True):
    """
    """

    atomic_symbols = []
    xyz_coordinates = []
    charge = 0
    title = ""

    with open(filename, "r") as file:
        for line_number, line in enumerate(file):
            if line_number == 0:
                num_atoms = int(line)
            elif line_number == 1:
                title = line
                if "charge=" in line:
                    charge = int(line.split("=")[1])
            else:
                atomic_symbol, x, y, z = line.split()
                atomic_symbols.append(atomic_symbol)
                xyz_coordinates.append([float(x), float(y), float(z)])

    atoms = [int_atom(atom) for atom in atomic_symbols]

    return atoms, charge, xyz_coordinates


def xyz2AC(atoms, xyz, charge, use_huckel=False, tr_previous_AC = [], N2collision = False):
    """

    atoms and coordinates to atom connectivity (AC)

    args:
        atoms - int atom types
        xyz - coordinates
        charge - molecule charge

    optional:
        use_huckel - Use Huckel method for atom connecitivty

    returns
        ac - atom connectivity matrix
        mol - rdkit molecule

    """

    if use_huckel:
        return xyz2AC_huckel(atoms, xyz, charge, N2collision = N2collision)
    else:
        return xyz2AC_vdW(atoms, xyz, tr_previous_AC, N2collision = N2collision)


def xyz2AC_vdW(atoms, xyz, tr_previous_AC = [], N2collision = False):

    # Get mol template
    mol = get_proto_mol(atoms)

    # Set coordinates
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, (xyz[i][0], xyz[i][1], xyz[i][2]))
    mol.AddConformer(conf)

    AC = get_AC(mol, covalent_factor = 1.40, tr_previous_AC = tr_previous_AC,
                N2collision=N2collision)

    return AC, mol


def get_AC(mol, covalent_factor=1.3, tr_previous_AC = [], N2collision=False):
    """

    Generate adjacent matrix from atoms and coordinates.

    AC is a (num_atoms, num_atoms) matrix with 1 being covalent bond and 0 is not


    covalent_factor - 1.3 is an arbitrary factor

    args:
        mol - rdkit molobj with 3D conformer

    optional
        covalent_factor - increase covalent bond length threshold with facto

    returns:
        AC - adjacent matrix

    """
    # Calculate distance matrix
    dMat = Chem.Get3DDistanceMatrix(mol)

    pt = Chem.GetPeriodicTable()
    num_atoms = mol.GetNumAtoms()
    AC = np.zeros((num_atoms, num_atoms), dtype=int)

    for i in range(num_atoms):
        a_i = mol.GetAtomWithIdx(i)
        Rcov_i = pt.GetRcovalent(a_i.GetAtomicNum()) * covalent_factor
        for j in range(i + 1, num_atoms):
            #this check is to not allow connections between the colliding N2 and the molecule
            if not N2collision or (i > num_atoms - 3) == (j > num_atoms - 3):
                a_j = mol.GetAtomWithIdx(j)
                Rcov_j = pt.GetRcovalent(a_j.GetAtomicNum()) * covalent_factor
                if (dMat[i, j] <= Rcov_i + Rcov_j):
                    AC[i, j] = 1
                    AC[j, i] = 1
                elif len(tr_previous_AC) > 0:
                    if tr_previous_AC[i, j] > 0 and dMat[i, j] <= (Rcov_i + Rcov_j) * 1.5:
                        AC[i, j] = 2
                        AC[j, i] = 2

                if AC[i, j] > 0:
                    # Check if total bonds for atom i is exceeded
                    if np.count_nonzero(AC[i, :]) > atomic_valence_electrons[a_i.GetAtomicNum()]:
                        # Find the bond with the highest distance and remove it
                        max_dist_j = np.argmax(dMat[i, :] * AC[i, :])
                        AC[i, max_dist_j] = 0
                        AC[max_dist_j, i] = 0

                        # Check if total bonds for atom j is exceeded
                    if np.count_nonzero(AC[j, :]) > atomic_valence_electrons[a_j.GetAtomicNum()]:
                        # Find the bond with the highest distance and remove it
                        max_dist_i = np.argmax(dMat[j, :] * AC[j, :])
                        AC[j, max_dist_i] = 0
                        AC[max_dist_i, j] = 0

    return AC


def xyz2AC_huckel(atomicNumList, xyz, charge, N2collision=False):
    """

    args
        atomicNumList - atom type list
        xyz - coordinates
        charge - molecule charge

    returns
        ac - atom connectivity
        mol - rdkit molecule

    """
    mol = get_proto_mol(atomicNumList)

    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i,(xyz[i][0],xyz[i][1],xyz[i][2]))
    mol.AddConformer(conf)

    num_atoms = len(atomicNumList)
    AC = np.zeros((num_atoms,num_atoms)).astype(int)

    mol_huckel = Chem.Mol(mol)
    mol_huckel.GetAtomWithIdx(0).SetFormalCharge(charge) #mol charge arbitrarily added to 1st atom    

    passed,result = rdEHTTools.RunMol(mol_huckel)
    opop = result.GetReducedOverlapPopulationMatrix()
    tri = np.zeros((num_atoms, num_atoms))
    tri[np.tril(np.ones((num_atoms, num_atoms), dtype=bool))] = opop #lower triangular to square matrix
    for i in range(num_atoms):
        for j in range(i+1,num_atoms):
            pair_pop = abs(tri[j,i])   
            if pair_pop >= 0.15: #arbitry cutoff for bond. May need adjustment
                AC[i,j] = 1
                AC[j,i] = 1

    return AC, mol


def chiral_stereo_check(mol):
    """
    Find and embed chiral information into the model based on the coordinates

    args:
        mol - rdkit molecule, with embeded conformer

    """
    Chem.SanitizeMol(mol)
    Chem.DetectBondStereochemistry(mol, -1)
    Chem.AssignStereochemistry(mol, flagPossibleStereoCenters=True, force=True)
    Chem.AssignAtomChiralTagsFromStructure(mol, -1)

    return


def xyz2mol(atoms, coordinates, charge=0, allow_charged_fragments=True,
            use_graph=True, use_huckel=False, embed_chiral=True,
            use_atom_maps=False, tr_previous_AC = [], N2collision = False):
    """
    Generate a rdkit molobj from atoms, coordinates and a total_charge.

    args:
        atoms - list of atom types (int)
        coordinates - 3xN Cartesian coordinates
        charge - total charge of the system (default: 0)

    optional:
        allow_charged_fragments - alternatively radicals are made
        use_graph - use graph (networkx)
        use_huckel - Use Huckel method for atom connectivity prediction
        embed_chiral - embed chiral information to the molecule

    returns:
        mols - list of rdkit molobjects

    """

    # Get atom connectivity (AC) matrix, list of atomic numbers, molecular charge,
    # and mol object with no connectivity information
    AC, mol = xyz2AC(atoms, coordinates, charge,
                     use_huckel=use_huckel, tr_previous_AC=tr_previous_AC,
                     N2collision=N2collision)

    # Convert AC to bond order matrix and add connectivity and charge info to
    # mol object
    new_mols, AC = AC2mol(mol, AC, atoms, charge,
                     allow_charged_fragments=allow_charged_fragments,
                     use_graph=use_graph,
                     use_atom_maps=use_atom_maps)

    # Check for stereocenters and chiral centers
    if embed_chiral:
        for new_mol in new_mols:
            chiral_stereo_check(new_mol)

    return new_mols, AC





if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(usage='%(prog)s [options] molecule.xyz')
    parser.add_argument('structure', metavar='structure', type=str)
    parser.add_argument('-s', '--sdf',
        action="store_true",
        help="Dump sdf file")
    parser.add_argument('--ignore-chiral',
        action="store_true",
        help="Ignore chiral centers")
    parser.add_argument('--no-charged-fragments',
        action="store_true",
        help="Allow radicals to be made")
    parser.add_argument('--no-graph',
        action="store_true",
        help="Run xyz2mol without networkx dependencies")

    # huckel uses extended Huckel bond orders to locate bonds (requires RDKit 2019.9.1 or later)
    # otherwise van der Waals radii are used
    parser.add_argument('--use-huckel',
        action="store_true",
        help="Use Huckel method for atom connectivity")
    parser.add_argument('-o', '--output-format',
        action="store",
        type=str,
        help="Output format [smiles,sdf] (default=sdf)")
    parser.add_argument('-c', '--charge',
        action="store",
        metavar="int",
        type=int,
        help="Total charge of the system")

    args = parser.parse_args()

    # read xyz file
    filename = args.structure

    # allow for charged fragments, alternatively radicals are made
    charged_fragments = not args.no_charged_fragments

    # quick is faster for large systems but requires networkx
    # if you don't want to install networkx set quick=False and
    # uncomment 'import networkx as nx' at the top of the file
    quick = not args.no_graph

    # chiral comment
    embed_chiral = not args.ignore_chiral

    # read atoms and coordinates. Try to find the charge
    atoms, charge, xyz_coordinates = read_xyz_file(filename)

    # huckel uses extended Huckel bond orders to locate bonds (requires RDKit 2019.9.1 or later)
    # otherwise van der Waals radii are used
    use_huckel = args.use_huckel

    # if explicit charge from args, set it
    if args.charge is not None:
        charge = int(args.charge)

    # Get the molobjs
    mols = xyz2mol(atoms, xyz_coordinates,
        charge=charge,
        use_graph=quick,
        allow_charged_fragments=charged_fragments,
        embed_chiral=embed_chiral,
        use_huckel=use_huckel)

    # Print output
    for mol in mols:
        if args.output_format == "sdf":
            txt = Chem.MolToMolBlock(mol)
            print(txt)

        else:
            # Canonical hack
            isomeric_smiles = not args.ignore_chiral
            smiles = Chem.MolToSmiles(mol, isomericSmiles=isomeric_smiles)
            m = Chem.MolFromSmiles(smiles)
            smiles = Chem.MolToSmiles(m, isomericSmiles=isomeric_smiles)
            print(smiles)
