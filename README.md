# `cebeconf`

```
             _                                __ 
            | |                              / _|
   ___  ___ | |__    ___   ___  ___   _ __  | |_ 
  / __|/ _ \| '_ \  / _ \ / __|/ _ \ | '_ \ |  _|
 | (__|  __/| |_) ||  __/| (__| (_) || | | || |  
  \___|\___||_.__/  \___| \___|\___/ |_| |_||_|
```

`cebeconf` package is a set of machine-learning models for predicting 1s-`c`ore `e`lectron `b`inding `e`nergies of `CONF` atoms in organic molecules (Ref-1). 

# Details of target-level 1s core-electron binding energies
- Models were trained on 12880 small organic molecules from the [bigQM7ω dataset](https://moldis-group.github.io/bigQM7w/) (Ref-2).
- Target property (1s core-electron binding energies) was calculated using the meta-GGA-DFT method strongly constrained and appropriately normed (`SCAN`) with a large, `Tight-full` numeric atom-centered orbital (NAO) basis set implemented in [FHI-aims](https://fhi-aims.org/).
- These calculations were performed using ωB97XD/def2TZVP geometries presented in the bigQM7ω dataset.
- For delta learning, the baseline energies were assigned based on Mulliken occupations. The data can be found in `Baseline_files`.
- Two example files (UFF-PBE : [ethane](https://github.com/moldis-group/cebeconf/blob/main/example_Mulliken_ethane_UFF_pbe_cc-pVDZ.txt) and [propane](https://github.com/moldis-group/cebeconf/blob/main/example_Mulliken_propane_UFF_pbe_cc-pVDZ.txt)) are also provided in home folder showing the output from Mulliken.out file from FHI-aims.

 # Details of training the ML models 
- To facilitate rapid application of the ML models, training was done using _baseline_ geometries of the bigQM7ω molecules determined with the universal force field (UFF). These geometries are also provided at [https://moldis-group.github.io/bigQM7w/](https://moldis-group.github.io/bigQM7w/)
- So, for new predictions, the ML models require geometries quickly determined with UFF.
- ML models were trained using the kernel-ridge-regression model using the atomic Coulomb matrix representation.
- For technical details, see Ref-1, and its Supporting Information. 

# Run `cebeconf` 

 - Install dependencies `numpy`, `pandas`

- Download and install the package
```
    git clone git@github.com:moldis-group/cebeconf.git
    pip3 install -e cebeconf
```
- Install from PyPI
```
   pip3 install cebeconf
```

 - Create an XYZ file at the UFF level (see below to learn about how to do this)

 - Run the ML model in `python3` (example in `cebeconf/test` folder)

 ```
from cebeconf import calc_be
  
calc_be('test.xyz','direct', 'ACM')
 ```

 - Suppose `test.xyz' contains the following geometry (which is the last molecule in bigQM7ω dataset)
```
18
bigQM7w_UFF_012883
C     1.03070  -0.07680   0.06770  
C     2.53800  -0.21440  -0.12550  
C     2.99750  -0.46340  -1.49170  
N     3.09380   0.90540  -0.90860  
C     4.47940   1.20090  -0.50870  
C     5.01760   2.53370  -1.00430  
C     4.47490   2.41010   0.41050  
H     0.59860  -1.07330   0.29480  
H     0.52630   0.33730  -0.83250  
H     0.83500   0.60170   0.92380  
H     3.17550  -0.57150   0.71420  
H     2.25180  -0.44020  -2.31440  
H     3.99580  -0.93590  -1.63370  
H     5.09800   0.43550   0.01500  
H     4.34280   2.85880  -1.82600  
H     6.09080   2.33310  -1.20820  
H     3.60210   3.09770   0.43410  
H     5.35240   2.60380   1.06330 
```

- Running the code generates the following output
```
...
 +--------------+
 | User inputs: |
 +--------------+
 Reading coordinates from: test.xyz
 Predicting 1s CEBEs using direct ML with the ACM descriptor 

 +--------------+
 | Prediction:  |
 +--------------+
    1 C      1.03070000     -0.07680000      0.06770000     290.81 eV
    2 C      2.53800000     -0.21440000     -0.12550000     291.83 eV
    3 C      2.99750000     -0.46340000     -1.49170000     291.90 eV
...
```

# How to calculate UFF-level geometry? 

Write down the [SMILES descriptor](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) of the molecule (example `c1ccccc1` for benzene) in a file. 

    echo 'c1ccccc1' > benzene.smi

Generate an initial geometry using [openbabel](http://openbabel.org/wiki/Main_Page). If you have obtained an initial geometry by other means, then you can skip the previous step.

    obabel -oxyz benzene.smi > benzene.xyz --gen3d

Relax tightly using UFF.

    obminimize -oxyz -sd -ff UFF -c 1e-8 benzene.xyz > benzene_UFF.xyz

:warning: We have used Open Babel 2.4.1 in our workflow.

:warning: We have tested the code to work with the following versions of Python and other modules. 
```
Python version: 3.8.20 
Numpy version: 1.21.6
PyTorch version: 2.4.1
SchNetPack version: 0.3
```

# References
[Ref-1] [_Chemical Space-Informed Machine Learning Models for
Rapid Predictions of X-ray Photoelectron Spectra of Organic Molecules_](https://doi.org/10.1088/2632-2153/ad871d)    
Susmita Tripathy, Surajit Das, Shweta Jindal, Raghunathan Ramakrishnan      
Mach. Learn.: Sci. Technol. 5 (2024) 045023.

[Ref-2] [_The Resolution-vs.-Accuracy Dilemma in Machine Learning Modeling of Electronic Excitation Spectra_](https://doi.org/10.1039/D1DD00031D)                  
Prakriti Kayastha, Sabyasachi Chakraborty, Raghunathan Ramakrishnan    
Digital Discov., 1 (2022) 689-702.    


