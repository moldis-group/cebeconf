# `cebeconf` :construction:

`cebeconf` is a set of machine-learning models of 1s-`c`ore `e`lectron `b`inding `e`nergies of `CONF` atoms in small organic molecules. 

# Details of target-level 1s core-electron binding energies
- Models were trained on 12880 small organic molecules from the [bigQM7ω dataset](https://moldis-group.github.io/bigQM7w/) (Ref-1).
- Target property, which is the 1s core-electron binding energies, were calculated using the meta-GGA-DFT method strongly constrained and appropriately normed (`SCAN`) with a `Tight-full` numeric atom-centered orbital (NAO) basis set implemented in the quantum mechanics software package, [FHI-aims](https://fhi-aims.org/).
- It is important to note that these calculations were performed using ωB97XD/def2TZVP geometries presented initially in the [bigQM7ω dataset](https://doi.org/10.1039/D1DD00031D), See [https://moldis-group.github.io/bigQM7w/](https://moldis-group.github.io/bigQM7w/).

 # Further details 
- To facilitate rapid application of the ML models, training was done using geometries of the bigQM7ω molecules [determined using the universal force field (UFF)](https://ndownloader.figshare.com/files/30478326) provided at [https://moldis-group.github.io/bigQM7w/](https://moldis-group.github.io/bigQM7w/)
- So, for new predictions, the ML models require geometries determined with UFF. 
- Additional technical details are summarized in an upcoming article, see Ref-3 below. 

# Install the package with `pip`

# How to quickly generate atomic coordinates with UFF?

Write down the [SMILES descriptor](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) of the molecule (example `c1ccccc1` for benzene) in a file. 

    echo 'c1ccccc1' > benzene.smi

Generate an initial geometry using [openbabel](http://openbabel.org/wiki/Main_Page). :information_desk_person: If you have obtained an initial geometry by other means, then you can skip the previous step.

    obabel -oxyz benzene.smi > benzene.xyz --gen3d

Relax tightly using UFF.

    obminimize -sd -ff UFF -c 1e-8 benzene.xyz > benzene_UFF.xyz

:warning: We have used Open Babel 2.4.1 in our workflow.

# Descriptor degeneration
(https://www.qmlcode.org/) (Ref-2).

# References
[Ref-1] [_The Resolution-vs.-Accuracy Dilemma in Machine Learning Modeling of Electronic Excitation Spectra_](https://doi.org/10.1039/D1DD00031D)                  
Prakriti Kayastha, Sabyasachi Chakraborty, Raghunathan Ramakrishnan    
Digital Discovery, 1 (2022) 689-702.    

[Ref-2] [_AS Christensen, FA Faber, B Huang, LA Bratholm, A Tkatchenko, KR Muller, OA von Lilienfeld (2017) "QML: A Python Toolkit for Quantum Machine Learning, https://github.com/qmlcode/qml"_](https://github.com/qmlcode/qml)  

[Ref-3][_Accurate Core-Electron Binding Energies using Machine Learning Models Trained on the Small Organic Molecules Chemical Space_](arxiv link)    
Susmita Tripathy, Shweta Jindal, Raghunathan Ramakrishnan      
To be posted in Arxiv. 
