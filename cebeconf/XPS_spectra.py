#Getting XPS spectra from core binding energies
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.stats as stats

#___________________PLOT BEs AS LINES________________________ #plot energies as lines of height(intensity) 1
def BE_lines(BE,color_element_i,ax_i): #C and BE are 1D arrays of length n
	y_new = 0.465
	for i in range(len(BE)):
		ax_i.axvline(x=BE[i], ymin=0.0, ymax=y_new, color=color_element_i[i], linestyle="--")
		if len(BE) != 1:
			ax_i.text(BE[i], y_new+0.15, str(i+1),ha="center")
		else:
			ax_i.text(BE[i], y_new-0.15, str(i+1),ha="center",va="bottom")
	return

#__________________CONVERT ENERGIES TO GAUSSIANS____________
def gaussian_func(BE,x,std_dev):
	gaussian = [0 for i in range(len(BE))]
	for i in range(len(BE)):
		gaussian[i] = stats.norm.pdf(x,BE[i],std_dev)
	return gaussian

#__________________SUM THE GAUSSIANS________________________
def make_spectra(x,gaussian,ax_i):
	XPS = [0 for i in range(len(gaussian[0]))]

	for i in range(len(XPS)):
		for j in range(len(gaussian)):
			XPS[i] += gaussian[j][i]
	ax_i.plot(x,XPS, color="black")
	return XPS

#__________________PLOTTING THE SPECTRA______________________
def plot_spectra(BE_element,color_element_i,ax_i,name_i):
	C = [i+1 for i in range(len(BE_element))]

	x = np.linspace(max(BE_element)+2,min(BE_element)-2,200)

	BE_lines(BE_element,color_element_i,ax_i)

	gaussian = gaussian_func(BE_element, x,0.7)
	for i in range(len(gaussian)):
		ax_i.plot(x, gaussian[i],color=color_element_i[i])#, label='$C_{}(XCH)$'.format({i+1}))'


	gaussian_KT = gaussian_func(BE_element,x,0.7)
	if len(C) > 1:
		XPS_KT = make_spectra(x, gaussian_KT,ax_i)

	ax_i.title.set_text(f"{name_i}")
	ax_i.set_xlabel(r"Binding energy (eV)")
	ax_i.set_ylabel(r"Intensity (arbritary units)")
	return

#___________________CALC_BE___________________________
def get_spectra(KRR_model, BE, BE_KT=0):
    if  KRR_model.lower() == 'direct':
        BE_KT = [0 for energy in BE]
    elif KRR_model.lower() == 'delta':
        if BE_KT == 0 or len(BE_KT) != len(BE):
            print("For Delta-ML, correct baseline energies were not provided as input, spectra will not be given as output. Add the delta-ML predictions to the baseline energies separately.")
            return
    BE = [energy1+energy2 for energy1,energy2 in zip(BE,BE_KT)]
    BE_atomwise = [[] for i in range(4)]

    for energy in BE:
        if energy > 250 and energy < 325:
            BE_atomwise[0].append(energy)
        elif energy > 350 and energy < 500:
            BE_atomwise[1].append(energy)
        elif energy > 500 and energy < 600:
            BE_atomwise[2].append(energy)
        elif energy > 600 and energy < 700:
            BE_atomwise[3].append(energy)
        else : print("ENERGY OUT OF RANGE!!",energy)
    fig, ax= plt.subplots(2,2)
    
    name = ["carbon", "nitrogen", "oxygen", "fluorine"]
    color_scale = [cm.Greys, cm.Blues, cm.Reds, cm.Greens]
    
    for i in range(4):
        if len(BE_atomwise[i]) == 0:
            ax[int(np.floor(i/2))][i%2].set_visible(False)
        else:
            plot_spectra(BE_atomwise[i],color_scale[i](np.linspace(0.25,0.75, len(BE_atomwise[i]))),ax[int(np.floor(i/2))][i%2],name[i])
    plt.tight_layout()
    plt.show()
    return
