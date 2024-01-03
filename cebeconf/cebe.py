import os
import numpy as np
import pandas as pd
import qml
from qml.representations import generate_atomic_coulomb_matrix
from setuptools import find_packages, setup
from pkg_resources import resource_filename

data_folder = resource_filename('cebeconf', 'data')

logo='''
             _                                __
            | |                              / _|
   ___  ___ | |__    ___   ___  ___   _ __  | |_
  / __|/ _ \| '_ \  / _ \ / __|/ _ \ | '_ \ |  _|
 | (__|  __/| |_) ||  __/| (__| (_) || | | || |
  \___|\___||_.__/  \___| \___|\___/ |_| |_||_|
'''


def calc_be(XYZfile):

    # Atomic numbers
    def atno(ele):
        if ele == 'H':
            Z = 1
        elif ele == 'C':
            Z = 6
        elif ele == 'N':
            Z = 7
        elif ele == 'O':
            Z = 8
        elif ele == 'F':
            Z = 9
        return Z

    # Main code
    print('')
    print(logo)
    print('')
    print(' This is an ML model for predicting 1s core binding energies calculated using the')
    print(' Delta-SCF approach with the metaGGA-DFT method, SCAN and a very large basis set.')
    print('')
    print(' Here are some standard values calculated with this DFT model')
    print('')
    print(' C in CH4, methane      290.94 eV')
    print(' C in CH3CH3, ethane    290.78 eV')
    print(' C in CH2CH2, ethylene  290.86 eV')
    print(' C in HCCH, acetylene   291.35 eV')
    print(' N in NH3               405.79 eV')
    print(' O in H2O               540.34 eV')
    print(' F in HF                694.95 eV')
    print('')

    # Read XYZfile
    mol_R=[]
    mol_Z=[]

    iline=0
    at_types=[]

    with open(XYZfile,'r') as f:

        print(' Reading geometry from', XYZfile)
        print('')

        for line in f:

            line=line.split()

            if iline == 0:
                N_at=int(line[0])

            if iline > 1:
                at_R=[float(line[1]),float(line[2]),float(line[3])]
                at_R=np.array(at_R)
                mol_R.append(at_R)

                ele=line[0]
                at_types.append(ele)
                at_Z=atno(ele)
                mol_Z.append(at_Z)

            iline=iline+1

    mol_Z = np.array(mol_Z)
    mol_R = np.array(mol_R)

    # Calculate descriptor for query molecule
    desc_q = generate_atomic_coulomb_matrix(mol_Z, mol_R, size=23, sorting='distance', central_cutoff=10.0, interaction_cutoff=10.0)

    # Load data
    X_train_C=np.load(os.path.join(data_folder, 'C_representation.npy'))
    X_train_N=np.load(os.path.join(data_folder, 'N_representation.npy'))
    X_train_O=np.load(os.path.join(data_folder, 'O_representation.npy'))
    X_train_F=np.load(os.path.join(data_folder, 'F_representation.npy'))

    df = pd.read_csv(os.path.join(data_folder, 'C_model_direct.csv'), header=None)
    model_C=np.array(df.iloc[:,0].values)
    df = pd.read_csv(os.path.join(data_folder, 'N_model_direct.csv'), header=None)
    model_N=np.array(df.iloc[:,0].values)
    df = pd.read_csv(os.path.join(data_folder, 'O_model_direct.csv'), header=None)
    model_O=np.array(df.iloc[:,0].values)
    df = pd.read_csv(os.path.join(data_folder, 'F_model_direct.csv'), header=None)
    model_F=np.array(df.iloc[:,0].values)

    sigma_C=5712.74014896
    sigma_N=9203.50735350
    sigma_O=12841.51702904
    sigma_F=87500.54782720

    # Gaussian and Laplacian kernels
    def kernel(option,sigma,dT, dQ):
        if option == 'L':
            dij=np.sum(np.abs(dT-dQ))
            val = np.exp(-dij / sigma)
        elif option == 'G':
            dij=np.sqrt(np.sum(np.abs(dT-dQ)**2))
            val = np.exp( -dij**2   / (2*sigma**2) )
        return val

    # Predict with KRR
    for i_at in range(N_at):

        avail=[6,7,8,9]

        if mol_Z[i_at] in avail:

            dQ=desc_q[i_at]
            dQ=np.array([dQ])

            if mol_Z[i_at] == 6:
                sigma=sigma_C
            elif mol_Z[i_at] == 7:
                sigma=sigma_N
            elif mol_Z[i_at] == 8:
                sigma=sigma_O
            elif mol_Z[i_at] == 9:
                sigma=sigma_F

            if mol_Z[i_at] == 6:

                Kpred=[]
                for i in range(len(X_train_C)):
                    dT=X_train_C[i]
                    Kiq=kernel('L',sigma,dT,dQ)
                    Kpred.append(Kiq)
                Epred=np.dot(Kpred,model_C)

            elif mol_Z[i_at] == 7:

                Kpred=[]
                for i in range(len(X_train_N)):
                    dT=X_train_N[i]
                    Kiq=kernel('L',sigma,dT,dQ)
                    Kpred.append(Kiq)
                Epred=np.dot(Kpred,model_N)

            elif mol_Z[i_at] == 8:

                Kpred=[]
                for i in range(len(X_train_O)):
                    dT=X_train_O[i]
                    Kiq=kernel('L',sigma,dT,dQ)
                    Kpred.append(Kiq)
                Epred=np.dot(Kpred,model_O)

            elif mol_Z[i_at] == 9:

                Kpred=[]
                for i in range(len(X_train_F)):
                    dT=X_train_F[i]
                    Kiq=kernel('L',sigma,dT,dQ)
                    Kpred.append(Kiq)
                Epred=np.dot(Kpred,model_F)

            print(f" atom: {i_at+1:4d} ({at_types[i_at]}), {Epred:6.3f} eV")

        else:

            print(f" atom: {i_at+1:4d} ({at_types[i_at]})")
