def rcut(mol_R, i_at, Z_i,dict_MaxN):
    '''
    Calculates cut-off radius adaptively for chopping large molecules with > 23 atoms

            Input:
                    mol_R(np.array, float): Atomic coordinates (N_at,3)
                    i_at (int): Index of the query atom
                    Z_i(int): Atomic number of atom-i to decide the number of neighbors
                    LargeSystem(logical): To activate an approximation for large systems

            Returns:
                    cutoff (float): cut-off radius
    '''
    import numpy as np

    Max_neigh=23

    if Z_i == 6:
        Max_neigh=dict_MaxN.get("MaxN_C")
    elif Z_i == 7:
        Max_neigh=dict_MaxN.get("MaxN_N")
    elif Z_i == 8:
        Max_neigh=dict_MaxN.get("MaxN_O")
    elif Z_i == 9:
        Max_neigh=dict_MaxN.get("MaxN_F")

    Ri=mol_R[i_at]
    cutoff=10.0

    N_at=len(mol_R)
    Nneigh=N_at

    while Nneigh > Max_neigh:

        Nneigh=0

        for j_at in range(N_at):

            Rj=mol_R[j_at]
            dRij=Ri-Rj
            Rij=np.sqrt(np.sum(dRij**2))

            if Rij < cutoff:

                Nneigh=Nneigh+1

        cutoff=cutoff-0.01

    return cutoff, Nneigh

def headers():
    '''
    Prepares header content for the output

            Input:

            Returns:
                    logo (string): cebeconf logo
                    header (string): header content
    '''

    logo='''
             _                                __
            | |                              / _|
   ___  ___ | |__    ___   ___  ___   _ __  | |_
  / __|/ _ \| '_ \  / _ \ / __|/ _ \ | '_ \ |  _|
 | (__|  __/| |_) ||  __/| (__| (_) || | | || |
  \___|\___||_.__/  \___| \___|\___/ |_| |_||_|
    '''

    header='''
 This is an ML model for predicting 1s core binding
 energies of CONF atoms. The model is trained on data
 calculated using Delta-SCF approach with the mGGA-DFT
 method, SCAN, at the basis set limit. 

 The model assumes that the input geometry is relaxed using
 universal force fields (UFF).

 If Delta-ML prediction is requested, the ML estimates
 the difference in SCAN Delta-SCF values from Koopmans
 estimation based on PBE/cc-pVDZ on UFF geometry.

 Some reference values determined with the SCAN Delta-SCF  method:

 C in CH4, methane      290.94 eV
 C in CH3CH3, ethane    290.78 eV
 C in CH2CH2, ethylene  290.86 eV
 C in HCCH, acetylene   291.35 eV
 N in NH3               405.79 eV
 O in H2O               540.34 eV
 F in HF                694.95 eV
    '''
    return logo, header

