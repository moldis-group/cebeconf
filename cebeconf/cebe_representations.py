import numpy as np

def damp(x, Rcut):
    val = 1.0 - 1.0/( 1.0 + np.exp(-(x-Rcut)))
    return val
  
#def LocalCM(mol_Z, mol_R, Max_at, Rcut):
def LocalCM(mol_Z, mol_R, Max_at):
    '''
    Calculates Local Coulomb Matrix for each atom
    '''
    indices=[]
    N_at=len(mol_Z)

    for i in range(N_at):

        distances = []
        for j in range(N_at):
            if i != j:
                distance = np.linalg.norm(mol_R[i] - mol_R[j])
                distances.append((j, distance))

        distances.sort(key=lambda x: x[1], reverse=False)

        indi=[]
        indi.append(i)
        for idx, dist in distances:
            indi.append(idx)

        indices.append(indi)

    desc=[]
    for k in range(N_at):

       #CMmat=np.zeros( [Max_at,Max_at] )
       #for i1,i in enumerate(indices[k]):
       #    rik = np.linalg.norm(mol_R[i] - mol_R[k])
       #    fik=damp(rik,Rcut)
       #    for j1,j in enumerate(indices[i]):
       #        rij = np.linalg.norm(mol_R[i] - mol_R[j])
       #        fij=damp(rij,Rcut)
       #        rjk = np.linalg.norm(mol_R[j] - mol_R[k])
       #        fjk=damp(rjk,Rcut)
       #        if i == j:
       #            val=0.5*mol_Z[i]**2.4
       #            CMmat[i1][i1]=val * fik**2
       #        else:
       #            val=mol_Z[i]*mol_Z[j]/rij
       #            CMmat[i1][j1]=val * fij * fjk * fik

        CMmat=np.zeros( [Max_at,Max_at] )
        for i1,i in enumerate(indices[k]):
            rik = np.linalg.norm(mol_R[i] - mol_R[k])
            for j1,j in enumerate(indices[i]):
                rij = np.linalg.norm(mol_R[i] - mol_R[j])
                if i == j:
                    val=0.5*mol_Z[i]**2.4
                    CMmat[i1][i1]=val
                else:
                    val=mol_Z[i]*mol_Z[j]/rij
                    CMmat[i1][j1]=val

        CM=[]
        for i in range(Max_at):
            for j in range(i,Max_at):
                CM.append(CMmat[i][j])
        desc.append(CM)

    desc=np.array(desc)

    return desc

def AtomicEnvt(mol_Z, mol_R, version=0.3):
    '''
    Calculates Local atomic environment for each atom
    '''
    import torch
    import schnetpack as spk
    from ase.io import read
    import os
    import warnings
    from pkg_resources import resource_filename

    if os.path.isfile('data.db') == True:
        os.remove('data.db')
    if os.path.isfile('data.pt') == True:
        os.remove('data.pt')
    if os.path.isfile('absurd_name.xyz') == True:
        os.remove('absurd_name.xyz')

    warnings.filterwarnings("ignore", message="The given NumPy array is not writable") #Numpy array [0.] is generated here, hence the warning
    warnings.filterwarnings("ignore", message="Use get_global_number_of_atoms()")  #Generated at for batch in test_loader
    
    atom_dic = {6:("C","model_C",[]), 7:("N","model_N",[]), 8:("O","model_O",[]), 9:("F","model_F",[]), 1:("H","DUMMY",[])}
    data_folder = resource_filename('cebeconf', 'models')


    #Dummy xyz as input for Schnet model
    f = open('absurd_name.xyz','w')
    print(len(mol_Z), file=f)
    print("eng=  0.000",end="", file=f)
    idx = 0
    desc = [[] for Z_i in mol_Z]
    for Z_i, R_i in zip(mol_Z,mol_R):
        if Z_i not in atom_dic:
            print("Atom Model not present.")
            exit
        else:
            print("\n"+atom_dic[Z_i][0], end="    ", file=f)
            for coord in R_i:
                print(coord, end="   ", file=f)
        atom_dic[Z_i][2].append(idx)
        idx += 1
    f.close()

    #Convert input file into .db file
    atoms = read('absurd_name.xyz', index=":")
    torch.save(atoms, "data.pt")
    data = torch.load("data.pt")
    property_ls =  [{"energy": np.array([mol.info["eng"]], dtype="float32")} for mol in data]
    
    dataset = spk.AtomsData("data.db", available_properties=["energy"])
    #dataset = spk.data.AtomsDataModule( 'data.db', batch_size=1, property_units={'energy':'eV'})
    dataset.add_systems(data, property_ls)
    test_loader = spk.AtomsLoader(dataset, batch_size=1)
    for item in atom_dic:
        if item != 1 and len(atom_dic[item][2]) > 0:
            model = torch.load(os.path.join(data_folder, atom_dic[item][1]), map_location=torch.device('cpu'))
            #model = torch.load(os.path.join(data_folder, 'model_C'), map_location=torch.device('cpu'))
            #print(f"{len(atom_dic[item][2])} {atom_dic[item][0]} atoms are present. {atom_dic[item][1]} will be used to generate representation ")
            for batch in test_loader:
                pred=model(batch)
                representation = batch["representation"]
                for idx in atom_dic[item][2]:
                    #print(representation[0][idx])
                    for number in representation[0][idx]:
                        #print(float(number.item()),end=" ")
                        desc[idx].append(number.item())
                    #print("")
                    #desc[idx] = torch.tensor.detach(representation[0][idx]).numpy()
        else:
            for idx in atom_dic[item][2]:
                desc[idx] = [0.0 for i in range(128)]
    #elif version == 3:
    #    print("Update this")
    desc = np.array(desc)
    #print(desc[0][0])
    os.remove('absurd_name.xyz')
    os.remove("data.db")
    os.remove("data.pt")
    return desc

