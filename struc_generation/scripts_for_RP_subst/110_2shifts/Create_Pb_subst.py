from ase.io import read,write
from ase.build.supercells import make_supercell
from ase import Atom, Atoms
import numpy as np

for (single_vac, double_vac) in [(1,0), (0,1)]: #(2, 0), (3, 0), (4, 0), (0, 2), (0, 3), (0, 4), (1,1), (2, 2)]:
    atoms = read('pure_struc.xyz')

    rand_lst_Pb = np.random.permutation(48)[:single_vac+double_vac]
    rand_lst_Pb_s = rand_lst_Pb[:single_vac]

    #Find I atoms with only one Pb neighbor
    Pb_dct = {}
    for i,at in enumerate(atoms):
        if at.symbol == "I":
            num_Pb_neighbors = 0
            j_Pb = None
            for j,at2 in enumerate(atoms):
                if at2.symbol == "Pb":
                    dist = atoms.get_distance(i, j, mic=True, vector=False)
                    if dist < 3.5:
                        j_Pb = j
                        num_Pb_neighbors += 1
            if num_Pb_neighbors == 1:
                if j_Pb not in Pb_dct.keys():
                    Pb_dct[j_Pb] = [i]
                else:
                    Pb_dct[j_Pb].append(i)

    tel = -1
    remove_I_lst = []
    remove_Cs_lst = []
    for i, at in enumerate(atoms):
        if at.symbol == 'Pb':
            tel += 1
            if tel in rand_lst_Pb:
                at.symbol = 'Cd'
                rand_int = np.random.randint(2)
                remove_I_lst.append(Pb_dct[i][rand_int])
                rand_lst = np.random.permutation(len(atoms))
                for k in rand_lst:
                    at3 = atoms[k]
                    if at3.symbol == 'Cs' and k not in remove_Cs_lst:
                        dist = atoms.get_distance(Pb_dct[i][rand_int], k, mic=True, vector=False)
                        dist_x1 = atoms.get_distance(Pb_dct[i][rand_int], i, mic=True, vector=True)[0]
                        dist_x2 = atoms.get_distance(Pb_dct[i][rand_int], k, mic=True, vector=True)[0]
                        if dist < 4.5 and np.sign(dist_x1) != np.sign(dist_x2):
                            remove_Cs_lst.append(k)
                            break
                if tel not in rand_lst_Pb_s:
                    remove_I_lst.append(Pb_dct[i][(rand_int +1) % 2])
                    rand_lst = np.random.permutation(len(atoms))
                    for k in rand_lst:
                        at3 = atoms[k]
                        if at3.symbol == 'Cs' and k not in remove_Cs_lst:
                            dist = atoms.get_distance(Pb_dct[i][(rand_int +1) % 2], k, mic=True, vector=False)
                            dist_x1 = atoms.get_distance(Pb_dct[i][(rand_int +1) % 2], i, mic=True, vector=True)[0]
                            dist_x2 = atoms.get_distance(Pb_dct[i][(rand_int +1) % 2], k, mic=True, vector=True)[0]
                            if dist < 4.5 and np.sign(dist_x1) != np.sign(dist_x2):
                                remove_Cs_lst.append(k)
                                break

    assert len(remove_I_lst) == len(remove_Cs_lst), "Number of I and Cs to remove is not equal"

    sym_lst = []
    pos_lst= []
    num_at_dct = {}
    for el in ["Cd", "Cs", "Pb", "I"]:
        count = 0
        for i, at in enumerate(atoms):
            if i not in remove_I_lst and i not in remove_Cs_lst:
                if at.symbol == el:
                    count += 1
                    sym_lst.append(at.symbol)
                    pos_lst.append(at.position.copy())
        num_at_dct[el] = count
    new_atoms = Atoms(symbols=sym_lst, positions=pos_lst, pbc=True, cell=atoms.cell)
    print("net charge: ", num_at_dct["Cs"] + 2*num_at_dct["Pb"] + 2*num_at_dct["Cd"] - num_at_dct["I"])

    if single_vac == 1 and double_vac == 0:
        write('atoms_110_Cd_Pb5co.xyz', new_atoms, format='extxyz')
    elif single_vac == 0 and double_vac == 1:
        write('atoms_110_Cd_Pb4co.xyz', new_atoms, format='extxyz')
    else:
        write('RP_110_n1_Cd_'+str(single_vac)+'_'+str(double_vac)+'_Ivac.xyz', new_atoms, format='extxyz')
        write('RP_110_n1_Cd_'+str(single_vac)+'_'+str(double_vac)+'_Ivac.cif', new_atoms, format='cif')

    for at in new_atoms:
        if at.symbol == "Cd":
            at.symbol = "Zn"
    if single_vac == 1 and double_vac == 0:
        write('atoms_110_Zn_Pb5co.xyz', new_atoms, format='extxyz')
    elif single_vac == 0 and double_vac == 1:
        write('atoms_110_Zn_Pb4co.xyz', new_atoms, format='extxyz')
    else:
        write('RP_110_n1_Zn_'+str(single_vac)+'_'+str(double_vac)+'_Ivac.xyz', new_atoms, format='extxyz')
        write('RP_110_n1_Zn_'+str(single_vac)+'_'+str(double_vac)+'_Ivac.cif', new_atoms, format='cif')

    for at in new_atoms:
        if at.symbol == "Zn":
            at.symbol = "Pb"
    if single_vac == 1 and double_vac == 0:
        write('atoms_110_Pb_Pb5co.xyz', new_atoms, format='extxyz')
    elif single_vac == 0 and double_vac == 1:
        write('atoms_110_Pb_Pb4co.xyz', new_atoms, format='extxyz')
    else:
        write('RP_110_n1_Pb_'+str(single_vac)+'_'+str(double_vac)+'_Ivac.xyz', new_atoms, format='extxyz')
        write('RP_110_n1_Pb_'+str(single_vac)+'_'+str(double_vac)+'_Ivac.cif', new_atoms, format='cif')

for num_Pb in [1]: #range(1,5):
    atoms = read('pure_struc.xyz')
    rand_lst_Pb = np.random.permutation(48)[:num_Pb]

    tel = -1
    for at in atoms:
        if at.symbol == 'Pb':
            tel += 1
            if tel in rand_lst_Pb:
                at.symbol = 'Cd'
    if num_Pb ==1:
        write('atoms_110_Cd_Pb6co.xyz', atoms, format='extxyz')
    else:
        write('RP_110_n1_Cd_'+str(num_Pb)+'.xyz', atoms, format='extxyz')
        write('RP_110_n1_Cd_'+str(num_Pb)+'.cif', atoms, format='cif')

    for at in atoms:
        if at.symbol == "Cd":
            at.symbol = "Zn"
    if num_Pb ==1:
        write('atoms_110_Zn_Pb6co.xyz', atoms, format='extxyz')
    else:
        write('RP_110_n1_Zn_'+str(num_Pb)+'.xyz', atoms, format='extxyz')
        write('RP_110_n1_Zn_'+str(num_Pb)+'.cif', atoms, format='cif')
