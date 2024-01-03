# `cebeconf` :construction:

```
             _                                __ 
            | |                              / _|
   ___  ___ | |__    ___   ___  ___   _ __  | |_ 
  / __|/ _ \| '_ \  / _ \ / __|/ _ \ | '_ \ |  _|
 | (__|  __/| |_) ||  __/| (__| (_) || | | || |  
  \___|\___||_.__/  \___| \___|\___/ |_| |_||_|
```

`cebeconf` package is a set of machine-learning models for predicting 1s-`c`ore `e`lectron `b`inding `e`nergies of `CONF` atoms in organic molecules. 

# Details of target-level 1s core-electron binding energies
- Models were trained on 12880 small organic molecules from the [bigQM7ω dataset](https://moldis-group.github.io/bigQM7w/) (Ref-1).
- Target property (1s core-electron binding energies) was calculated using the meta-GGA-DFT method strongly constrained and appropriately normed (`SCAN`) with a large, `Tight-full` numeric atom-centered orbital (NAO) basis set implemented in [FHI-aims](https://fhi-aims.org/).
- These calculations were performed using ωB97XD/def2TZVP geometries presented in the bigQM7ω dataset.

 # Details of training the ML models 🤖
- To facilitate rapid application of the ML models, training was done using _baseline_ geometries of the bigQM7ω molecules determined with the universal force field (UFF). These geometries are also provided at [https://moldis-group.github.io/bigQM7w/](https://moldis-group.github.io/bigQM7w/)
- So, for new predictions, the ML models require geometries quickly determined with UFF.
- ML models were trained using the kernel-ridge-regression model using the atomic Coulomb matrix representation (Ref-2).
- Additional technical details are summarized in an upcoming article (Ref-3). 

# Run `cebeconf` ✨

 - Install dependencies `numpy`, `qml`, `pandas`

- Download and install the package
```
    mkdir Test_cebeconf
    cd Test_cebeconf
    git clone git@github.com:moldis-group/cebeconf.git
    pip3 install -e cebeconf
```
 - Create an XYZ file at the UFF level (see below to learn about how to do this)

 - Run the ML model in `python3` (example in `cebeconf/test` folder)

 ```
from cebeconf import calc_be
  
calc_be('test.xyz')
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
 C      1.03070000     -0.07680000      0.06770000     290.62 eV, 0.24 seconds
 C      2.53800000     -0.21440000     -0.12550000     291.22 eV, 0.26 seconds
 C      2.99750000     -0.46340000     -1.49170000     291.04 eV, 0.26 seconds
 N      3.09380000      0.90540000     -0.90860000     404.81 eV, 0.11 seconds
...
```

# How to calculate UFF-level geometry? 🤔

Write down the [SMILES descriptor](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) of the molecule (example `c1ccccc1` for benzene) in a file. 

    echo 'c1ccccc1' > benzene.smi

Generate an initial geometry using [openbabel](http://openbabel.org/wiki/Main_Page). :information_desk_person: If you have obtained an initial geometry by other means, then you can skip the previous step.

    obabel -oxyz benzene.smi > benzene.xyz --gen3d

Relax tightly using UFF.

    obminimize -sd -ff UFF -c 1e-8 benzene.xyz > benzene_UFF.xyz

:warning: We have used Open Babel 2.4.1 in our workflow.

# References
[Ref-1] [_The Resolution-vs.-Accuracy Dilemma in Machine Learning Modeling of Electronic Excitation Spectra_](https://doi.org/10.1039/D1DD00031D)                  
Prakriti Kayastha, Sabyasachi Chakraborty, Raghunathan Ramakrishnan    
Digital Discovery, 1 (2022) 689-702.    

[Ref-2] [_AS Christensen, FA Faber, B Huang, LA Bratholm, A Tkatchenko, KR Muller, OA von Lilienfeld (2017) "QML: A Python Toolkit for Quantum Machine Learning, https://github.com/qmlcode/qml"_](https://github.com/qmlcode/qml)  

[Ref-3][_Accurate Core-Electron Binding Energies using Machine Learning Models Trained on the Small Organic Molecules Chemical Space_](arxiv link)    
Susmita Tripathy, Shweta Jindal, Raghunathan Ramakrishnan      
To be posted in Arxiv. 
