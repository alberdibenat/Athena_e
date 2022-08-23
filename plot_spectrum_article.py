import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#photon_flux = pd.read_csv("SR_paper_final.dc0",skiprows=[0,1,2,3,4,5,6,7,8,9], header=None, delim_whitespace=True)
photon_flux = pd.read_csv("Simulation_updated_article-8.txt",skiprows=[0,1], header=None, delim_whitespace=True)
##photon_flux.columns=['Energy','Flux','PL', 'PC', 'PL45','1-PL']
photon_flux.columns=['Energy','Flux','PL', 'PC', 'PL45']

#coherent_photon_flux = pd.read_csv("Simulation_updated_article-5.txt",skiprows=[0,1,2,3,4,5,6,7,8,9], header=None, delim_whitespace=True)
coherent_photon_flux = pd.read_csv("Simulation_CRS_athenae-10.txt",skiprows=[0,1], header=None, delim_whitespace=True)
coherent_photon_flux.columns=['Energy','Flux','PL', 'PC', 'PL45']
h = 6.62e-34
e = 1.6e-19
c = 3e8


energy = np.array(photon_flux['Energy'].astype(float).tolist()) #photon energy in eV
flux = np.array(photon_flux['Flux'].astype(float).tolist()) #photons/s/0.1%
coherent_energy = np.array(coherent_photon_flux['Energy'].astype(float).tolist()) #photon energy in eV
coherent_flux = np.array(coherent_photon_flux['Flux'].astype(float).tolist()) #photons/s/0.1%
wavelength = h*c/(energy*e)
coherent_wavelength = h*c/(coherent_energy*e)



wavelength_plot = np.flip(wavelength)
coherent_wavelength_plot = np.flip(coherent_wavelength)
flux_plot = np.flip(flux)
coherent_flux_plot = np.flip(coherent_flux)



flux_per_20nm = np.array([])
flux_per_50nm = np.array([])
flux_per_105nm = np.array([])
for i in range(len(energy)):
    wavelength_value = wavelength[i]
    BW = wavelength_value*0.001
    flux_per_20nm = np.append(flux_per_20nm,flux[i]*20e-9/BW)
    flux_per_50nm = np.append(flux_per_50nm,flux[i]*50e-9/BW)
    flux_per_105nm = np.append(flux_per_105nm,flux[i]*105e-9/BW)

coherent_flux_per_20nm = np.array([])
coherent_flux_per_50nm = np.array([])
coherent_flux_per_105nm = np.array([])
for i in range(len(coherent_energy)):
    coherent_wavelength_value = coherent_wavelength[i]
    coherent_BW = coherent_wavelength_value*0.001
    coherent_flux_per_20nm = np.append(coherent_flux_per_20nm,coherent_flux[i]*20e-9/coherent_BW)
    coherent_flux_per_50nm = np.append(coherent_flux_per_50nm,coherent_flux[i]*50e-9/coherent_BW)
    coherent_flux_per_105nm = np.append(coherent_flux_per_105nm,coherent_flux[i]*105e-9/coherent_BW)


colors = ['blue','red','green','darkorange','darkgrey']
plt.rcParams['axes.labelsize'] = 24
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['font.weight']='bold'
plt.rcParams['font.size'] = 20
plt.rcParams['legend.handlelength'] = 2 #Added by me
plt.rcParams['legend.markerscale'] = 2 #Added by me



fig2,ax21 = plt.subplots(figsize=(11,8))
#ax21.set_xlabel(r'$\lambda$ [nm]')
ax21.set_xlabel(r'$\lambda$ [nm]')
ax21.set_ylabel('Photon flux @ detector (ph/s)')
#ax21.set_title(r'SR in bending dipole')
ax21.plot(np.flip(c*h/energy/e)*1e9,np.flip(flux_per_20nm),linewidth='3',color=colors[0],label='ISR cSTART, BW = 20 nm')#'SR 910mT 20nm BW')
ax21.plot(np.flip(c*h/energy/e)*1e9,np.flip(flux_per_50nm),linewidth='3',color=colors[3],label='ISR cSTART, BW = 50 nm')#'SR 910mT 50nm BW')
ax21.plot(np.flip(c*h/energy/e)*1e9,np.flip(flux_per_105nm),linewidth='3',color=colors[2],label='ISR cSTART, BW = 105 nm')#'SR 910mT 100nm BW')
ax21.plot(np.flip(c*h/coherent_energy/e)*1e9,np.flip(coherent_flux_per_20nm),color=colors[0],linewidth='3',linestyle='--',label=r'CSR $\mathrm{Athena_e}$, BW = 20 nm')#'SR 910mT 20nm BW')
ax21.plot(np.flip(c*h/coherent_energy/e)*1e9,np.flip(coherent_flux_per_50nm),color=colors[3],linewidth='3',linestyle='--',label=r'CSR $\mathrm{Athena_e}$, BW = 50 nm')#'SR 910mT 50nm BW')
ax21.plot(np.flip(c*h/coherent_energy/e)*1e9,np.flip(coherent_flux_per_105nm),color=colors[2],linewidth='3',linestyle='--',label=r'CSR $\mathrm{Athena_e}$, BW = 105 nm')#'SR 910mT 100nm BW')
#ax21.plot(coherent_energy*e/h,coherent_flux,label='CSR 910mT 100cm')
ax21.set_ylim(0.01,5e11)
#ax21.set_xlim(1e-3*e/h,30*e/h) #for frequency
ax21.set_xlim(c/(30*e/h)*1e9,c/(1e-3*e/h)*1e9) #for wavelength
#ax21.axvspan(380.0e12, 750.0e12, alpha=0.3, color='grey') #for frequency
ax21.axvspan(c/750.0e12*1e9,c/380.0e12*1e9, alpha=0.3, color='grey') #for wavelength
#ax21.plot(wavelength_plot*1e9,flux_plot)
#ax21.text(0.4, 0.8, 'Visible light', rotation=-90, va='center')
ax21.text(550, 100, 'Visible light',
         rotation= -90,
         horizontalalignment='center',
         verticalalignment='center',
         multialignment='center')
plt.grid()
plt.xscale('log')
plt.yscale('log')
fig2.tight_layout()
plt.legend(loc='upper right')
fig2.savefig('SR_dipole_article.jpg',dpi=1200)
plt.show()
