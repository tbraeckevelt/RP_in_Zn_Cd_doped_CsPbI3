from ase.io import read,write
import numpy as np
from ase import Atom, Atoms


for num_Cd in [1,3,6]:
    atoms = read("pure_struc.xyz")

    #get random list of numbers between 0 and 96 with length num_Cd
    rand_int_lst = np.random.permutation(96)[:num_Cd]

    tel = -1
    ind_Cs_lst = []
    Cs_lst = []
    for i, at in enumerate(atoms):
        if at.symbol == 'Cs':
            Cs_lst.append(i)
            tel += 1
            if tel in rand_int_lst:
                ind_Cs_lst.append(i)

    #find Cs atoms that are within 6 Angstroms of the randomly selected Cs atoms
    Cs_pairs_lst = []
    remove_lst = []
    for ind_Cs in ind_Cs_lst:
        for i in Cs_lst:
            dist = atoms.get_distance(ind_Cs, i, mic=True, vector=False)
            if dist < 6 and i not in ind_Cs_lst and i not in remove_lst:
                dist_vec = atoms.get_distance(ind_Cs, i, mic=True, vector=True)
                if np.abs(dist_vec[2]) > 1 and np.abs(dist_vec[1]) + np.abs(dist_vec[0]) > 1:
                    Cs_pairs_lst.append((ind_Cs, i))
                    remove_lst.append(i)
                    break
    
    print(Cs_pairs_lst)
    print(remove_lst)
    print(ind_Cs_lst)

    #create structure with Cs atoms removed and Cd atoms added 
    sym_lst = []
    pos_lst= []
    for el in ["Cs", "Pb", "I"]:
        for i, at in enumerate(atoms):
            if i not in ind_Cs_lst and i not in remove_lst:
                if at.symbol == el:
                    sym_lst.append(at.symbol)
                    pos_lst.append(at.position.copy())
    for (i,j) in Cs_pairs_lst:
        dist_vec = atoms.get_distance(i, j, mic=True, vector=True)
        pos_lst.append(atoms[i].position.copy() + dist_vec/2)
        sym_lst.append("Cd")

    new_atoms = Atoms(symbols=sym_lst, positions=pos_lst, pbc=True, cell=atoms.cell)
    write("Cd_subst_"+str(num_Cd)+".xyz", new_atoms)