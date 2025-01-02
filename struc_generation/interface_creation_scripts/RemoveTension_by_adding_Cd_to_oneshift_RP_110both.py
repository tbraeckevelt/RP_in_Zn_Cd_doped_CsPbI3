from ase.io import read,write
from ase import Atom, Atoms
import numpy as np


num_adap = 32
half_RP_x = 20   #needed later
atoms = read("Perovskite_RP_110both_big_varying_oneshift.xyz")

#create list with all the indices of a speficied element
ind_dct = {"Cs":[], "Pb":[], "I":[]}
for i, at in enumerate(atoms):
    ind_dct[at.symbol].append(i)

#create list of all I atoms with only 1 Pb neighbor

I1_lst = []
for i in ind_dct["I"]:
    num_Pb = 0
    for j in ind_dct["Pb"]:
        dist = atoms.get_distance(i, j, mic=True)
        if dist < 3.5:
            num_Pb += 1
    assert num_Pb == 1 or num_Pb == 2, "I atom has " +str(num_Pb)+ " Pb neighbors"
    if num_Pb == 1:
        I1_lst.append(i)
        print(i)
print(I1_lst)

#adaption option 3: change Pb at the edges to Cd and change their I atoms coordination to 4 instead of 6
Pb_edge_lst = []
I_edge_lst = []
for i in ind_dct["Pb"]:
    I1_neig = []
    for j in I1_lst:
        dist = atoms.get_distance(i, j, mic=True)
        if dist < 3.5:
            I1_neig.append(j)
    assert len(I1_neig) < 5, "Pb atom has " +str(len(I1_neig))+ " I1 neighbors"
    if len(I1_neig) > 3:
        Pb_edge_lst.append(i)
        for j in I1_neig:
            tel = 0
            for j2 in I1_neig:
                if j != j2:
                    dist = atoms.get_distance(j, j2, mic=True, vector=True)
                    if np.abs(dist[0]) > 4.0:
                        tel += 1
            if tel != 2:
                I_edge_lst.append(j)

rand_lst_Pb = np.random.permutation(len(Pb_edge_lst))[:num_adap]
remove_lst = []
tel = -1
for i in Pb_edge_lst:
    tel += 1
    if tel in rand_lst_Pb:
        atoms[i].symbol = "Cd"
        for j in I_edge_lst:
            dist = atoms.get_distance(i, j, mic=True)
            if dist < 3.5:
                for k in ind_dct["Cs"]:
                    dist = atoms.get_distance(j, k, mic=True)
                    if dist < 3.5 and j not in remove_lst and k not in remove_lst:
                        remove_lst.append(j)
                        remove_lst.append(k)
                        break

sym_lst = []
pos_lst = []
num_at_dct = {}
for el in ["Cs", "Pb", "I", "Cd"]:
    count = 0
    for at in atoms:
        if at.symbol == el and at.index not in remove_lst:
            count +=1
            sym_lst.append(at.number)
            pos_lst.append(at.position.copy())
    num_at_dct[el] = count
print(num_at_dct)   
print("net charge option 3: ", num_at_dct["Cs"] + 2*num_at_dct["Pb"] + 2*num_at_dct["Cd"] - num_at_dct["I"])

new_atoms = Atoms(symbols=sym_lst, positions=pos_lst, pbc=True, cell=atoms.cell)

#outname = "testingstruc_option3_"+str(num_adap)+".cif"
#write(outname, new_atoms, format="cif")
outname = "Perovskite_RP_110both_big_varying_oneshift_Cd"+str(num_adap)+".xyz"
write(outname, new_atoms, format="extxyz")


#The other option you can find below, we did not apply these options in the end because they did not work well.
#Option1: The Pb octahedra at the edge had to stretch to much which lead to breaking of the Pb-I bonds and distortion of the structure
#Option2: The Cd atoms at the interface could pull out the Pb atoms at the edge which lead to distortion of the structure

"""
#adaption option 1: merge closeby I atoms and Cs atoms
#Create list of all I atoms which have a second Pb neighbor at a distance between 3.5 and 5

I2_lst = []
for i in I1_lst:
    for j in ind_dct["Pb"]:
        dist = atoms.get_distance(i, j, mic=True)
        if dist > 3.5 and dist < 5:
            I2_lst.append(i)
            print(i)
            break
print(I2_lst)

#Make pairs of I atoms of I2_lst who are neighbors
pair_I_lst = []
remove_lst = []
for tel, i in enumerate(I2_lst):
    for j in I2_lst[tel:]:
        if i != j:
            dist = atoms.get_distance(i, j, mic=True)
            if dist < 3.5:
                pair_I_lst.append((i,j))
                print(i, j)
                remove_lst.append(i)
                remove_lst.append(j)
                break

#Find for each pair of I atoms a pair Cs atoms closeby
pair_Cs_lst = []
for i, j in pair_I_lst:
    for k in ind_dct["Cs"]:
        dist = atoms.get_distance(i, k, mic=True, vector=True)
        if (dist[0]**2 + dist[2]**2) < 0.1 and np.abs(dist[1]) < 3.5 and np.sign(atoms[i].position[0] - half_RP_x)*dist[1] < 0: 
            i1 = k
        elif np.abs(dist[1]) < 0.01 and (dist[0]**2 + dist[2]**2) < 3.5**2 and np.sign(atoms[i].position[0] - half_RP_x) * dist[0] > 0:
            i1 = k
        dist = atoms.get_distance(j, k, mic=True, vector=True)
        if (dist[0]**2 + dist[2]**2) < 0.1 and np.abs(dist[1]) < 3.5 and np.sign(atoms[i].position[0] - half_RP_x) * dist[1] < 0:
            i2 = k
        elif np.abs(dist[1]) < 0.01 and (dist[0]**2 + dist[2]**2) < 3.5**2 and np.sign(atoms[i].position[0] - half_RP_x) * dist[0] > 0:
            i2 = k
    pair_Cs_lst.append((i1, i2))
    print(i1, i2)
    remove_lst.append(i1)
    remove_lst.append(i2)

pos_lst = []
sym_lst = []
for i, at in enumerate(atoms):
    if at.symbol == "Cs":
        if i not in remove_lst:
            pos_lst.append(at.position.copy())
            sym_lst.append(at.symbol)
for i, j in pair_Cs_lst:
    dist = atoms.get_distance(i, j, mic=True, vector=True)
    pos_lst.append(atoms[i].position.copy() + dist/2)
    sym_lst.append("Cs")
for i, at in enumerate(atoms):
    if at.symbol == "Pb":
        pos_lst.append(at.position.copy())
        sym_lst.append(at.symbol)
for i, at in enumerate(atoms):
    if at.symbol == "I":
        if i not in remove_lst:
            pos_lst.append(at.position.copy())
            sym_lst.append(at.symbol)
for i, j in pair_I_lst:
    dist = atoms.get_distance(i, j, mic=True, vector=True)
    pos_lst.append(atoms[i].position.copy() + dist/2)
    sym_lst.append("I")

new_atoms = Atoms(symbols=sym_lst, positions=pos_lst, pbc=True, cell=atoms.cell)

outname = "testingstruc_option1.cif"
write(outname, new_atoms, format="cif")
outname = "ForMLP_Per_RP_110both_shift_option1.xyz"
write(outname, new_atoms, format="extxyz")

#create list with all the indices of a speficied element
new_ind_dct = {"Cs":[], "Pb":[], "I":[]}
for i, at in enumerate(new_atoms):
    new_ind_dct[at.symbol].append(i)
print("net charge option 1: ", len(new_ind_dct["Cs"]) + 2*len(new_ind_dct["Pb"]) - len(new_ind_dct["I"]))

#adaption option 2: remove two closeby Cs atoms with one Cd atoms at the avearge position
pos_lst = []
sym_lst = []
for i, at in enumerate(atoms):
    if at.symbol == "Cs":
        if i not in remove_lst:
            pos_lst.append(at.position.copy())
            sym_lst.append(at.symbol)
for i, j in pair_Cs_lst:
    dist = atoms.get_distance(i, j, mic=True, vector=True)
    pos_lst.append(atoms[i].position.copy() + dist/2)
    sym_lst.append("Cd")
for i, at in enumerate(atoms):
    if at.symbol == "Pb":
        pos_lst.append(at.position.copy())
        sym_lst.append(at.symbol)
for i, at in enumerate(atoms):
    if at.symbol == "I":
        pos_lst.append(at.position.copy())
        sym_lst.append(at.symbol)


new_atoms = Atoms(symbols=sym_lst, positions=pos_lst, pbc=True, cell=atoms.cell)

outname = "testingstruc_option2.cif"
write(outname, new_atoms, format="cif")
outname = "ForMLP_Per_RP_110both_shift_option2.xyz"
write(outname, new_atoms, format="extxyz")

#create list with all the indices of a speficied element
new_ind_dct = {"Cs":[], "Pb":[], "I":[], "Cd":[]}
for i, at in enumerate(new_atoms):
    new_ind_dct[at.symbol].append(i)
print("net charge option 2: ", len(new_ind_dct["Cs"]) + 2*(len(new_ind_dct["Pb"]) + len(new_ind_dct["Cd"])) - len(new_ind_dct["I"]))
"""