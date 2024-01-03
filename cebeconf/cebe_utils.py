def rcut(mol_R, i_at):
   '''
    Calculates cut-off radius adaptively for chopping large molecules with > 23 atoms

            Input:
                    mol_R(np.array, float): Atomic coordinates (N_at,3)
                    i_at (int): Index of the query atom

            Returns:
                    cutoff (float): cut-off radius
    '''
    import numpy as np

    Ri=mol_R[i_at]
    cutoff=10.0

    N_at=len(mol_R)
    Nneigh=N_at

    while Nneigh > 23:
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
 method, SCAN, and a very large basis set.

 Some reference values determined with this DFT method:

 C in CH4, methane      290.94 eV
 C in CH3CH3, ethane    290.78 eV
 C in CH2CH2, ethylene  290.86 eV
 C in HCCH, acetylene   291.35 eV
 N in NH3               405.79 eV
 O in H2O               540.34 eV
 F in HF                694.95 eV
    '''
    return logo, header

