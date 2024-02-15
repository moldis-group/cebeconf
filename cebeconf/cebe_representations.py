import numpy as np

def damp(x, Rcut):
    val = 1.0 - 1.0/( 1.0 + np.exp(-(x-Rcut)))
    return val
  
def LocalCM(mol_Z, mol_R, Max_at, Rcut):
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

        CMmat=np.zeros( [Max_at,Max_at] )
        for i in indices[k]:
            rik = np.linalg.norm(mol_R[i] - mol_R[k])
            fik=damp(rik,Rcut)
            for j in indices[i]:
                rij = np.linalg.norm(mol_R[i] - mol_R[j])
                fij=damp(rij,Rcut)
                rjk = np.linalg.norm(mol_R[j] - mol_R[k])
                fjk=damp(rjk,Rcut)
                if i == j:
                    val=0.5*mol_Z[i]**2.4
                    CMmat[i][i]=val*fik**2
                else:
                    val=mol_Z[i]*mol_Z[j]/rij
                    CMmat[i][j]=val*fij*fjk*fjk

        CM=[]
        for i in range(Max_at):
            for j in range(i,Max_at):
                CM.append(CMmat[i][j])
        desc.append(CM)

    desc=np.array(desc)

    return desc

