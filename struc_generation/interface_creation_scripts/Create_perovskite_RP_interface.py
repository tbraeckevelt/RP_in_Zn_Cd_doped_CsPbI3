from ase.io import read, write
from ase import Atom, Atoms
import numpy as np
from ase.build import niggli_reduce
from ase.build.supercells import make_supercell

atoms = Atoms()
atoms_dup = Atoms()

l = 3.2 # Pb-I distance
total_length =  14
depth = 6
RP_dct = [(1,9),(1,3)] #each instance is an RP perovskite layer with (thickness n, end (for even instances) or start (for odd instances) position)
RP_type = "110" #plane of the CsI layer of the RP phase
inter_type = "110" #plane of the CsI layer of the interface between the perovskite and the RP phase
shift =  np.array([1, 1, 1]) * l #Shift between perovskite layers, 
                                 #Only positive values are allowed, otherwise noshift list and/or cell definition must be adapted

#CsPbI3 building block
atoms_fu = Atoms(symbols=["Cs", "Cs", "Cs", "Cs", "Cs", "Cs", "Cs", "Cs", "Pb",  "I", "I", "I", "I", "I", "I"], 
                    positions=[np.array([ l, l, l]), 
                            np.array([-l, l, l]), 
                            np.array([ l,-l, l]), 
                            np.array([-l,-l, l]), 
                            np.array([ l, l,-l]), 
                            np.array([-l, l,-l]), 
                            np.array([ l,-l,-l]), 
                            np.array([-l,-l,-l]), 
                            np.array([0.0, 0.0, 0.0]), 
                            np.array([  l, 0.0, 0.0]), 
                            np.array([ -l, 0.0, 0.0]), 
                            np.array([0.0,   l, 0.0]), 
                            np.array([0.0,  -l, 0.0]), 
                            np.array([0.0, 0.0,   l]), 
                            np.array([0.0, 0.0,  -l])])

noshift_lst = []
for wl, (n_p, co) in enumerate(RP_dct):
    prev_co = RP_dct[(wl-1)%len(RP_dct)][1]
    if co > prev_co:
        noshift_lst.append([])
        for i in range(total_length):
            if i < prev_co:
                noshift_lst[-1].append(True)
            elif i < co:
                if RP_type == "110" and inter_type == "110" and i == co -1:
                    noshift_lst[-1].append(False)
                else:
                    noshift_lst[-1].append(None)
            else:
                noshift_lst[-1].append(False)
        if RP_type == "110":
            noshift_lst.append([])
            for i in range(total_length):
                if i < prev_co:
                    noshift_lst[-1].append(True)
                elif i < co:
                    if inter_type == "110" and i == prev_co :
                        noshift_lst[-1].append(True)
                    else:
                        noshift_lst[-1].append(None)
                else:
                    noshift_lst[-1].append(False)
    for j in range(n_p):
        noshift_lst.append([])
        for i in range(total_length):
            if i < co:
                noshift_lst[-1].append(True)
            else:
                noshift_lst[-1].append(False)
    if RP_type == "110":
        noshift_lst.append([])
        for i in range(total_length):
            if i < co:
                noshift_lst[-1].append(True)
            else:
                noshift_lst[-1].append(False)

cell_np= np.array([[2*l*total_length, 0.0, 0.0] + shift, [0.0, 2*l*len(noshift_lst), 0.0], [0.0, 0.0, 2*l]])
if inter_type == "110":
    cell_np[2,0] =  -2*l
if RP_type == "110":
    cell_np[2,1] =  -2*l

for j in range(len(noshift_lst)):
    for i in range(total_length):
        print(j,i,noshift_lst[j][i])
        if noshift_lst[j][i] is not None:
            for at in atoms_fu.copy():
                if noshift_lst[j][i]:
                    at.position += np.array([l*2*i, l*2*j, 0.0])
                else:
                    at.position += np.array([l*2*i, l*2*j, 0.0]) + shift
                atoms_dup.append(at)

#Remove duplicate atoms
for at in atoms_dup:
    dir_pos = np.dot(at.position, np.linalg.inv(cell_np))
    dir_pos -= np.floor(dir_pos+0.0001)
    at.position = np.dot(dir_pos, cell_np)
    flag = True
    for at2 in atoms:
        dist = np.linalg.norm(at.position - at2.position)
        if dist < 0.1:
            flag = False
            #assert at.symbol == at2.symbol
            if at.symbol != at2.symbol:
                print("Different symbol at same position("+str(at.position)+"):" + at.symbol +" and " + at2.symbol)
                at2.symbol = "I"
    if flag:
        atoms.append(at)

at_numb_lst = []
pos_lst = []
num_at_dct = {}
for el in ["Cs", "Pb", "I"]:
    count = 0
    for at in atoms:
        if at.symbol == el:
            count +=1
            at_numb_lst.append(at.number)
            pos_lst.append(at.position.copy())
    num_at_dct[el] = count
print(num_at_dct)   
print("net charge: ", num_at_dct["Cs"] + 2*num_at_dct["Pb"] - num_at_dct["I"])
if num_at_dct["Cs"] + 2*num_at_dct["Pb"] - num_at_dct["I"] != 0:
    print("Warning: net charge is not zero, some atoms are overlapping that should not overlap, probably due to the setted shift")

new_atoms = Atoms(numbers=at_numb_lst, positions=pos_lst, pbc=True, cell=cell_np)

sup_atoms = make_supercell(new_atoms, np.diag([1,1,depth]))
niggli_reduce(sup_atoms)

outname = "teststruc.cif"
write(outname, sup_atoms, format="cif")
#outname = "Perovskite_RP_110both_big_varying.xyz"
#write(outname, sup_atoms, format="extxyz")


