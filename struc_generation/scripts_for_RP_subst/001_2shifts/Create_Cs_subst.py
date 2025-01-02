from ase.io import read,write
import numpy as np
from ase import Atom, Atoms

atoms = read("pure_struc.xyz")

#get random integer between 0 and 96
rand_int = np.random.randint(96)

tel = -1
ind_Cs = 0
for i, at in enumerate(atoms):
    if at.symbol == 'Cs':
        tel += 1
        if tel == rand_int:
            ind_Cs = i

#find all Cs atoms that are within 6.5 Angstroms of the randomly selected Cs atom
Cs_pairs_lst = []
pair_ind_lst = [-1, -1, -1]
for i, at in enumerate(atoms):
    if at.symbol == 'Cs':
        dist = atoms.get_distance(ind_Cs, i, mic=True, vector=False)
        if dist < 6.5 and i != ind_Cs:
            Cs_pairs_lst.append(i)
            dist_vec = atoms.get_distance(ind_Cs, i, mic=True, vector=True)
            if dist_vec[2]< -1:
                pair_ind_lst[0] = i
            elif dist_vec[2]> 1:
                pair_ind_lst[1] = i
            else:
                pair_ind_lst[2] = i

assert len(Cs_pairs_lst) == 9, "did not find enough neighbours, found "+str(len(Cs_pairs_lst))+" instead of 9"

#create structure with 2 Cs atoms removed and a Zn or Cd atom (in two different ways) for each pair of Cs atoms
for j, pair_ind in enumerate(pair_ind_lst):
    sym_lst = []
    pos_lst= []
    for el in ["Cs", "Pb", "I"]:
        for i, at in enumerate(atoms):
            if i != ind_Cs and i != pair_ind:
                if at.symbol == el:
                    sym_lst.append(at.symbol)
                    pos_lst.append(at.position.copy())

    pos_lst.append(atoms[ind_Cs].position.copy())
    sym_lst.append("Zn")
    new_atoms = Atoms(symbols=sym_lst, positions=pos_lst, pbc=True, cell=atoms.cell)
    write("Zn_subst_"+str(j)+"_Cspos.xyz", new_atoms)
    sym_lst[-1] = "Cd"
    new_atoms = Atoms(symbols=sym_lst, positions=pos_lst, pbc=True, cell=atoms.cell)
    write("Cd_subst_"+str(j)+"_Cspos.xyz", new_atoms)

    dist_vec = atoms.get_distance(ind_Cs, pair_ind, mic=True, vector=True)
    pos_lst[-1] = atoms[ind_Cs].position.copy() + dist_vec/2
    new_atoms = Atoms(symbols=sym_lst, positions=pos_lst, pbc=True, cell=atoms.cell)
    write("Cd_subst_"+str(j)+"_ave.xyz", new_atoms)
    sym_lst[-1] = "Zn"
    new_atoms = Atoms(symbols=sym_lst, positions=pos_lst, pbc=True, cell=atoms.cell)
    write("Zn_subst_"+str(j)+"_ave.xyz", new_atoms)