from ase.io import read,write
import numpy as np
from ase import Atom, Atoms


for num_Cd in [1,3,6]:

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

    #get random list of numbers between 0 and 96 with length num_Cd
    rand_int_lst = np.random.permutation(48)[:num_Cd]
    ind_Cs_lst = []
    for i, ind_Cs in enumerate(Cs_lst_2):
        if i in rand_int_lst:
            ind_Cs_lst.append(ind_Cs)

    #find Cs atoms that are within 6 Angstroms of the randomly selected Cs atoms
    Cs_pairs_lst = []
    remove_lst = []
    for ind_Cs in ind_Cs_lst:
        for i in Cs_lst_2:
            dist = atoms.get_distance(ind_Cs, i, mic=True, vector=False)
            if dist < 6 and i not in ind_Cs_lst and i not in remove_lst:
                dist_vec = atoms.get_distance(ind_Cs, i, mic=True, vector=True)
                if np.abs(dist_vec[1]) > 1:
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
