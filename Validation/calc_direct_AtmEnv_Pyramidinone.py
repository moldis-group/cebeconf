#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np 
import os 
from cebeconf import calc_be

folders = ['1H_2pyramidinone','1H_4pyramidinone','3H_2pyramidinone','3H_4pyramidinone']
loop = 4

ie = []

for folder in folders:
    path = './pyramidinone/' + folder
#   os.chdir(path)

    folder_name = path +'/'+folder+ '_direct.dat'
    with open(folder_name, 'r') as file:
        comment = file.readline()
        dft_data = file.readlines() 
    file.close()

    count = 0 

    for p in range(loop):
        for q in range(loop):
            for r in range(loop):
                molecule_name = path +'/'+ folder + '_' + str(p) + '_' + str(q) + '_'                     + str(r) + '.xyz'

                out_acm = calc_be(molecule_name, 'direct', 'AtmEnv')

                with open(molecule_name, 'r') as file:
                    no_of_atoms = file.readline()
                    molecule_name = file.readline()
                    lines = file.readlines() 
                file.close()

                ie_1 = []

                for i in range(np.size(lines)):
                    ie_2 = []
                    if lines[i].split()[0] != 'H': 
                        a, x, y, z = lines[i].split()
                        dft = dft_data[count].split()[2] 

                        ie_2.append(str(i+1))
                        ie_2.append(a)
                        ie_2.append(x)
                        ie_2.append(y)
                        ie_2.append(z)
                        ie_2.append(out_acm[i])   
                        ie_2.append(dft)

                        ie_1.append(ie_2)
                        count = count + 1
                ie.append(ie_1)

                # ie [Molecule] [Atom] [Number,Symbole,x,y,z,DFT,IE]
                #print(ie) 


# In[ ]:


count = 0 

path = './pyramidinone' 
os.chdir(path)

with open('pyramidinone.txt','w') as file:
    file.write("Atom     x        y         z         IE(DFT)      IE(ML)\n")
    for f in range(4):
        for p in range(loop):
            for q in range(loop):
                for r in range(loop):
                    molecule_name = folders[f] + '_' + str(p) + '_' + str(q) + '_'                         + str(r) + '.xyz'
                    file.write(molecule_name[:-4] + '\n')

                    for i in range(np.shape(ie[count])[0]):
                        n = ie[count][i][0]
                        a = ie[count][i][1]
                        x = float(ie[count][i][2])
                        y = float(ie[count][i][3])
                        z = float(ie[count][i][4]) 
                        dft = float(ie[count][i][5])
                        ml = float(ie[count][i][6])
                        dIE = dft-ml

                        file.write("{} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:12.5f} {:10.5f}                             \n".format(a, x, y, z, dft, ml, dIE)
                            )
                    count = count + 1 


# In[ ]:


# Standard Deviation Calculation 
c_count = 0
n_count = 0
o_count = 0
f_count = 0

c_list = []
n_list = []
o_list = []
f_list = []

count = 0 

for f in range(4):
    for p in range(loop):
            for q in range(loop):
                for r in range(loop):
                    for i in range(np.shape(ie[count])[0]):

                        atom = ie[count][i][1] 
                        dft = float(ie[count][i][5])
                        ml = float(ie[count][i][6])


                        if atom == 'C':
                            c_count += 1 
                            delta = abs(dft-ml)
                            c_list.append(delta)
                        if atom == 'N':
                            n_count += 1 
                            delta = abs(dft-ml)
                            n_list.append(delta)
                        if atom == 'O':
                            o_count += 1 
                            delta = abs(dft-ml)
                            o_list.append(delta)
                        if atom == 'F':
                            f_count += 1 
                            delta = abs(dft-ml)
                            f_list.append(delta)
                    count = count + 1 


if c_count != 0:
    print("#C = {}   MAE = {:9.4f} SD = {:8.4f}".format(c_count,np.mean(c_list),np.std(c_list)))
if n_count != 0:
    print("#N = {}   MAE = {:9.4f} SD = {:8.4f}".format(n_count,np.mean(n_list),np.std(n_list)))
if o_count != 0:
    print("#O = {}   MAE = {:9.4f} SD = {:8.4f}".format(o_count,np.mean(o_list),np.std(o_list)))
if f_count != 0:
    print("#F = {}   MAE = {:9.4f} SD = {:8.4f}".format(f_count,np.mean(f_list),np.std(f_list)))


# In[ ]:




