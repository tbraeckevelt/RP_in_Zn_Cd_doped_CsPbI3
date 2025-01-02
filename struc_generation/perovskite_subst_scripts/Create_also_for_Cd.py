from ase.io import read,write


for (num_int_Cs, num_int_Pb) in [(1, 0), (2, 0), (3, 0), (4, 0), (5, 0), (6, 0), (7, 0), (8, 0), 
                                 (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), 
                                 (1, 1), (2, 2), (3, 3), (4, 4), (3, 5), (5, 3), (3, 4), (4, 3)]:

    atoms = read('Alpha_Zn_subst_'+str(num_int_Cs)+'_'+str(num_int_Pb)+'.xyz')
    for at in atoms:
        if at.symbol == 'Zn':
            at.symbol = 'Cd'
    write('Alpha_Cd_subst_'+str(num_int_Cs)+'_'+str(num_int_Pb)+'.xyz', atoms, format='extxyz')
