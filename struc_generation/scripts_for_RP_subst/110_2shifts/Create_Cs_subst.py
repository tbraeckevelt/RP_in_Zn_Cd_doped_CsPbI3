from ase.io import read,write
import numpy as np
from ase import Atom, Atoms

atoms = read("pure_struc.xyz")

I_lst_1 = []  #I atoms in perovskite layer
I_lst_2 = []  #I atoms in between layers

for i, at in enumerate(atoms):
    if at.symbol == 'I':
        Pb_neigh = 0
        for j,at2 in enumerate(atoms):
            if at2.symbol == 'Pb':
                dist = atoms.get_distance(i, j, mic=True, vector=False)
                if dist < 3.5:
                    dist_vec = atoms.get_distance(i, j, mic=True, vector=True)
                    if np.abs(dist_vec[2]) < 3:         #Do not add the I atoms that are at the same X and Y position as the Pb atoms
                        Pb_neigh +=1
        assert Pb_neigh <= 2, "I atom "+str(i)+" has "+str(Pb_neigh)+" Pb neighbours"
        if Pb_neigh == 2:
            I_lst_1.append(i)
        elif Pb_neigh == 1:
            I_lst_2.append(i)

assert len(I_lst_1) == 48, "did not find enough I atoms in perovskite layer, found "+str(len(I_lst_1))+" instead of 48"
assert len(I_lst_2) == 96, "did not find enough I atoms in between layers, found "+str(len(I_lst_2))+" instead of 96"

Cs_lst_1 = []  #Cs atoms in perovskite layer
Cs_lst_2 = []  #Cs atoms in between layers

for i, at in enumerate(atoms):
    if at.symbol == "Cs":
        flag = True
        for j,at2 in enumerate(atoms):
            if at2.symbol == 'I':
                dist = atoms.get_distance(i, j, mic=True, vector=False)
                if dist < 4.4:
                    if j in I_lst_1:
                        Cs_lst_1.append(i)
                        flag = False
                        break
        if flag:
            Cs_lst_2.append(i)

assert len(Cs_lst_1) == 48, "did not find enough Cs atoms in perovskite layer, found "+str(len(Cs_lst_1))+" instead of 48"
assert len(Cs_lst_2) == 48, "did not find enough Cs atoms in between layers, found "+str(len(Cs_lst_2))+" instead of 48"


#get random integer between 0 and 96
rand_int = np.random.randint(48)
rand_ind1 = Cs_lst_1[rand_int]
pair_lst1 = [-1,-1,-1,-1]

for i in Cs_lst_1:
    if i != rand_ind1:
        dist = atoms.get_distance(rand_ind1, i, mic=True, vector=False)
        if dist < 6.5:
            dist_vec = atoms.get_distance(rand_ind1, i, mic=True, vector=True)
            if np.abs(dist_vec[0]) > 1:
                pair_lst1[0] = i
            else:
                pair_lst1[1] = i

for i in Cs_lst_2:
    dist = atoms.get_distance(rand_ind1, i, mic=True, vector=False)
    if dist < 6.7:
        rand_ind2 = i
        dist_vec = atoms.get_distance(rand_ind1, i, mic=True, vector=True)
        if np.abs(dist_vec[1]) > 1:
            pair_lst1[2] = i
        else:
            pair_lst1[3] = i

pair_lst2 = [-1,-1]
for i in Cs_lst_2:
    if i != rand_ind2:
        dist = atoms.get_distance(rand_ind2, i, mic=True, vector=False)
        if dist < 6.5:
            dist_vec = atoms.get_distance(rand_ind2, i, mic=True, vector=True)
            if np.abs(dist_vec[1]) > 1:
                pair_lst2[0] = i
            else:
                pair_lst2[1] = i


#create structure with 2 Cs atoms removed and a Zn or Cd atom (in different ways) for each pair of Cs atoms
for j, pair_ind in enumerate(pair_lst1):
    sym_lst = []
    pos_lst= []
    for el in ["Cs", "Pb", "I"]:
        for i, at in enumerate(atoms):
            if i != rand_ind1 and i != pair_ind:
                if at.symbol == el:
                    sym_lst.append(at.symbol)
                    pos_lst.append(at.position.copy())

    dist_vec = atoms.get_distance(rand_ind1, pair_ind, mic=True, vector=True)
    pos_lst.append(atoms[rand_ind1].position.copy() + dist_vec/2)
    sym_lst.append("Zn")
    new_atoms = Atoms(symbols=sym_lst, positions=pos_lst, pbc=True, cell=atoms.cell)
    write("Zn_subst_"+str(j)+"_ave.xyz", new_atoms)
    sym_lst[-1] = "Cd"
    new_atoms = Atoms(symbols=sym_lst, positions=pos_lst, pbc=True, cell=atoms.cell)
    write("Cd_subst_"+str(j)+"_ave.xyz", new_atoms)
        
for j, pair_ind in enumerate(pair_lst2):
    sym_lst = []
    pos_lst= []
    for el in ["Cs", "Pb", "I"]:
        for i, at in enumerate(atoms):
            if i != rand_ind2 and i != pair_ind:
                if at.symbol == el:
                    sym_lst.append(at.symbol)
                    pos_lst.append(at.position.copy())

    dist_vec = atoms.get_distance(rand_ind2, pair_ind, mic=True, vector=True)
    pos_lst.append(atoms[rand_ind2].position.copy() + dist_vec/2)
    sym_lst.append("Zn")
    new_atoms = Atoms(symbols=sym_lst, positions=pos_lst, pbc=True, cell=atoms.cell)
    write("Zn_subst_"+str(j+4)+"_ave.xyz", new_atoms)
    sym_lst[-1] = "Cd"
    new_atoms = Atoms(symbols=sym_lst, positions=pos_lst, pbc=True, cell=atoms.cell)
    write("Cd_subst_"+str(j+4)+"_ave.xyz", new_atoms)
