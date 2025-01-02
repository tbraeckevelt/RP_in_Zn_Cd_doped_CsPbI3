from ase.io import read,write
from ase.build.supercells import make_supercell
from ase import Atom, Atoms
import numpy as np


atoms = read('Alpha.vasp')

for (num_int_Cs, num_int_Pb) in [(1, 0), (2, 0), (3, 0), (4, 0), (5, 0), (6, 0), (7, 0), (8, 0), 
                                 (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), 
                                 (1, 1), (2, 2), (3, 3), (4, 4), (3, 5), (5, 3), (3, 4), (4, 3)]:
    #Make supercell
    atoms_sup = make_supercell(atoms, 4 * np.identity(3))

    rand_lst_Cs = np.random.permutation(64)[:num_int_Cs]
    rand_lst_Pb = np.random.permutation(64)[:num_int_Pb]

    #Add interstitials
    i = -1
    j = -1
    el_types = []
    for at in atoms_sup:
        el_types.append(at.symbol)
        if at.symbol == 'Cs':
            i += 1
            if i in rand_lst_Cs:
                #Add Zn interstitial at Cs site
                at_Zn = Atom('Zn', position=at.position.copy() - np.array([1,1,1]))
                atoms_sup.append(at_Zn)

                #Shift position Cs and change to I
                at.position += np.array([1,1,1])
                at.symbol = 'I'
        if at.symbol == 'Pb':
            j += 1
            if j in rand_lst_Pb:
                at.symbol = 'Zn' 
    for i in range(num_int_Cs):
        el_types.append('Zn')

    at_numb_lst = []
    pos_lst = []
    for el in ["Zn", "Cs", "Pb", "I"]:
        for at in atoms_sup:
            if at.symbol == el:
                at_numb_lst.append(at.number)
                pos_lst.append(at.position.copy())
                    
    new_atoms = Atoms(numbers=at_numb_lst, positions=pos_lst, pbc=True, cell=atoms_sup.cell)

    write('Alpha_Zn_subst_'+str(num_int_Cs)+'_'+str(num_int_Pb)+'.xyz', new_atoms, format='extxyz')
        
