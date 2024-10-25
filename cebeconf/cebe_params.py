def atno(at):
    '''
    Returns atom number

            Input:
                    at (str): Atomic symbol

            Returns:
                    atno (int): Atom number 
    '''
    at_num_dict ={"H":  1,  "He":  2,  "Li":  3,  "Be":  4,   "B":  5,\
                  "C":  6,   "N":  7,   "O":  8,   "F":  9,  "Ne": 10,\
                 "Na": 11,  "Mg": 12,  "Al": 13,  "Si": 14,   "P": 15,\
                  "S": 16,  "Cl": 17,  "Ar": 18,   "K": 19,  "Ca": 20,\
                 "Sc": 21,  "Ti": 22,   "V": 23,  "Cr": 24,  "Mn": 25,\
                 "Fe": 26,  "Co": 27,  "Ni": 28,  "Cu": 29,  "Zn": 30,\
                 "Ga": 31,  "Ge": 32,  "As": 33,  "Se": 34,  "Br": 35,\
                 "Kr": 36,   "I": 53}

    if at in at_num_dict:
        atno = at_num_dict[at]
        return(atno)
    else:
        print('Error: Unsupported element encountered'+str(at))
