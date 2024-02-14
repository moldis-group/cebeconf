import torch
import schnetpack as spk
from ase.io import read
from numpy import array
from copy import deepcopy
import os
import warnings

warnings.filterwarnings("ignore", message="The given NumPy array is not writable") #Numpy array [0.] is generated here, hence the warning
warnings.filterwarnings("ignore", message="Use get_global_number_of_atoms()")  #Generated at for batch in test_loader

dic = {"C":([],"model_C"), "N":([],"model_N"), "O":([],"model_O"), "F":([],"model_F"), "H":([],"DUMMY")}

def get_representation(XYZfile):
    f = open(XYZfile,'r')
    Lines = f.readlines()
    orig_line = Lines[1]
    Lines[1] = "eng=  0.000\n"  #safeguard for empty name or molecule name
    
    #Each element: [idx], "model_element"
    for idx,line in enumerate(Lines[2:]):
        line=line.split()
        dic[line[0]][0].append(idx)
    
    f.close()
    f = open(XYZfile,'w')
    f.writelines(Lines)
    f.close()
    
    #Convert input file into .db file
    atoms = read(XYZfile, index=":")
    torch.save(atoms, "data.pt")
    data = torch.load("data.pt")
    property_ls =  [{"energy": array([mol.info["eng"]], dtype="float32")} for mol in data]
    dataset = spk.AtomsData("data.db", available_properties=["energy"])
    dataset.add_systems(data, property_ls)
    test_loader = spk.AtomsLoader(dataset, batch_size=1)

    for item in dic:
        if item != "H" and len(dic[item][0]) > 0:
            print(f"{len(dic[item][0])} {item} atoms are present, {dic[item][-1]} will be used to generate representation file, '{item}_rep.txt'.")
            model = torch.load("models/"+dic[item][-1], map_location=torch.device('cpu'))
            output_file = open(item+"_rep.txt", "w")
            for batch in test_loader:
                pred=model(batch)
                representation = batch["representation"]
                for idx in dic[item][0]:
                    for number in representation[0][idx][:-2]:
                        output_file.write(f"{number},")
                    output_file.write(f"{representation[0][idx][-1]}\n")
            output_file.close()
    
    Lines[1] = orig_line
    f = open(XYZfile,"w")
    f.writelines(Lines)
    f.close()

    os.remove("data.db")
    os.remove("data.pt")
    return
