from ase.io import read,write
import numpy as np
from ase import Atom, Atoms


for num_Pb in [1,3,6]:
    for coor in [6,5,4]:
        rand_lst_Pb = np.random.permutation(48)[:num_Pb]

        atoms = read("pure_struc.xyz")

        I_lst = []  #I atoms with one neighbour
        Pb_lst = [] 
        Cs_lst = []

        tel = -1
        for i, at in enumerate(atoms):
            if at.symbol == 'I':
                Pb_neigh = 0
                for j,at2 in enumerate(atoms):
                    if at2.symbol == 'Pb':
                        dist = atoms.get_distance(i, j, mic=True, vector=False)
                        if dist < 3.5:
                            Pb_neigh +=1
                assert Pb_neigh == 1 or Pb_neigh == 2, "I atom "+str(i)+" has "+str(Pb_neigh)+" Pb neighbours"
                if Pb_neigh == 1:
                    I_lst.append(i)
            elif at.symbol == 'Pb':
                tel += 1
                if tel in rand_lst_Pb:
                    Pb_lst.append(i)
            elif at.symbol == 'Cs':
                Cs_lst.append(i)

        assert len(I_lst) == 96, "did not find enough I atoms in between layers, found "+str(len(I_lst))+" instead of 96"

        remove_lst = []
        for i, at in enumerate(atoms):
            if i in Pb_lst:
                at.symbol = 'Cd'
                if coor != 6:
                    telj = 0
                    for j in I_lst:
                        dist = atoms.get_distance(i, j, mic=True, vector=False)
                        if dist < 3.5:
                            remove_lst.append(j)
                            telj+=1
                            telk = 0
                            for k in Cs_lst:
                                dist = atoms.get_distance(j, k, mic=True, vector=False)
                                if dist < 5 and k not in remove_lst:
                                    remove_lst.append(k)
                                    telk+=1
                                    break
                            assert telk == 1, "I atom "+str(j)+" has removed "+str(telk)+" Cs neighbours"
                            if coor == 5:
                                break
                    assert (coor==4 and telj==2) or (coor==5 and telj==1),"Pb atom "+str(i)+" has removed "+str(telj)+" I neighbours for coordination "+str(coor)
                
        sym_lst = []
        pos_lst= []
        for el in ["Cs", "Pb", "I", "Cd"]:
            for i, at in enumerate(atoms):
                if i not in remove_lst:
                    if at.symbol == el:
                        sym_lst.append(at.symbol)
                        pos_lst.append(at.position.copy())


        new_atoms = Atoms(symbols=sym_lst, positions=pos_lst, pbc=True, cell=atoms.cell)
        write("Cd_subst_co"+str(coor)+"_"+str(num_Pb)+".xyz", new_atoms)

        """
        for i, el in enumerate(sym_lst):
            if el == "Zn":
                sym_lst[i] = "Cd"        
        new_atoms = Atoms(symbols=sym_lst, positions=pos_lst, pbc=True, cell=atoms.cell)
        write("Cd_subst_"+str(num_Pb)+".xyz", new_atoms)
        """