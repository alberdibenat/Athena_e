import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from numpy.random import random
import csv
from datetime import datetime
from scipy.optimize import curve_fit
from uncertainties import ufloat

def sinc(x):
    return np.sin(x)/x

def sinc2(x):
    return (np.sin(x)/x)**2

def interference_pattern(position,visibility,phi):
    #Geometry of experimental layout
    a = 0.0005
    d = 0.001
    f = 0.1
    #Central wavelength (green light)
    wavelength = 550e-9
    x1 = np.pi*a/(wavelength*f)*(position)
    x2 = 2*np.pi*d/(wavelength*f)*(position+phi)
    return sinc2(x1)*(1+visibility*np.cos(x2))/(1+visibility) #Normalized to 1


bandwidth_values = [15.0,30.0,50.0,75.0,105.0,140.0]
photon_numbers = [20000,100000,1000000] 
V_values = [0.100855,0.230338,0.472799,0.692774,0.960037]
Recovered_V_n_BW = []
Uncertainty_V_n_BW = []
for j in range(len(V_values)):
    Recovered_V_n_BW.append([])
    Uncertainty_V_n_BW.append([])
    filename = 'Raw_data_2/Photons_V_'+str('%02d'%(V_values[j]*100))+'_raw_data.txt'   #Write the projection in x into a txt file, this will simulate an experimental readout of a screen with 2um pixels and 1e5 incident photons
    photon_dataframe = pd.read_csv(filename, delim_whitespace = True)
    photon_dataframe.columns=['15','30','50','75','105','140']
    nm15 = np.array(photon_dataframe['15'].astype(float).tolist())
    nm30 = np.array(photon_dataframe['30'].astype(float).tolist())
    nm50 = np.array(photon_dataframe['50'].astype(float).tolist())
    nm75 = np.array(photon_dataframe['75'].astype(float).tolist())
    nm105 = np.array(photon_dataframe['105'].astype(float).tolist())
    nm140 = np.array(photon_dataframe['140'].astype(float).tolist())
    
    for n in range(len(photon_numbers)):
        Recovered_V_n_BW[j].append([])
        Uncertainty_V_n_BW[j].append([])
        measurements = 50
        bin_edges = []
        value = []
        rms_spread = []
        V15_carrier = []
        V30_carrier = []
        V50_carrier = []
        V75_carrier = []
        V105_carrier = []
        V140_carrier = []
        for i in range(measurements):
            nm15_new = nm15[i*photon_numbers[n]:(i+1)*photon_numbers[n]]
            nm30_new = nm30[i*photon_numbers[n]:(i+1)*photon_numbers[n]]
            nm50_new = nm50[i*photon_numbers[n]:(i+1)*photon_numbers[n]]
            nm75_new = nm75[i*photon_numbers[n]:(i+1)*photon_numbers[n]]
            nm105_new = nm105[i*photon_numbers[n]:(i+1)*photon_numbers[n]]
            nm140_new = nm140[i*photon_numbers[n]:(i+1)*photon_numbers[n]]

            counts15, xedges = np.histogram(nm15_new,bins=100,range=[-0.0002, 0.0002]) #4um pixel size asumed in x (range -0.0002m to 0.0002m divided in 100 pixels)
            counts30, xedges = np.histogram(nm30_new,bins=100,range=[-0.0002, 0.0002]) #4um pixel size asumed in x (range -0.0002m to 0.0002m divided in 100 pixels)
            counts50, xedges = np.histogram(nm50_new,bins=100,range=[-0.0002, 0.0002]) #4um pixel size asumed in x (range -0.0002m to 0.0002m divided in 100 pixels)
            counts75, xedges = np.histogram(nm75_new,bins=100,range=[-0.0002, 0.0002]) #4um pixel size asumed in x (range -0.0002m to 0.0002m divided in 100 pixels)
            counts105, xedges = np.histogram(nm105_new,bins=100,range=[-0.0002, 0.0002]) #4um pixel size asumed in x (range -0.0002m to 0.0002m divided in 100 pixels)
            counts140, xedges = np.histogram(nm140_new,bins=100,range=[-0.0002, 0.0002]) #4um pixel size asumed in x (range -0.0002m to 0.0002m divided in 100 pixels)
        
            max15 = np.amax(counts15)
            max30 = np.amax(counts30)
            max50 = np.amax(counts50)
            max75 = np.amax(counts75)
            max105 = np.amax(counts105)
            max140 = np.amax(counts140)

            counts15 = counts15/max15
            counts30 = counts30/max30
            counts50 = counts50/max50
            counts75 = counts75/max75
            counts105 = counts105/max105
            counts140 = counts140/max140

            param_bounds = ([0,-np.pi],[1,np.pi])
            x_points = [(xedges[l+1]+xedges[l])/2.0 for l in range(len(xedges)-1)] 
            popt15,pcov15 = curve_fit(interference_pattern,x_points,counts15,p0=(0.7,0.0),bounds=param_bounds)
            popt30,pcov30 = curve_fit(interference_pattern,x_points,counts30,p0=(0.7,0.0),bounds=param_bounds)
            popt50,pcov50 = curve_fit(interference_pattern,x_points,counts50,p0=(0.7,0.0),bounds=param_bounds)
            popt75,pcov75 = curve_fit(interference_pattern,x_points,counts75,p0=(0.7,0.0),bounds=param_bounds)
            popt105,pcov105 = curve_fit(interference_pattern,x_points,counts105,p0=(0.7,0.0),bounds=param_bounds)
            popt140,pcov140 = curve_fit(interference_pattern,x_points,counts140,p0=(0.7,0.0),bounds=param_bounds)
    
            V15_carrier.append(popt15)
            V30_carrier.append(popt30)
            V50_carrier.append(popt50)
            V75_carrier.append(popt75)
            V105_carrier.append(popt105)
            V140_carrier.append(popt140)
    
        print(len(nm15_new))

        mean_value_V15 = np.sum(V15_carrier)/measurements
        mean_value_V30 = np.sum(V30_carrier)/measurements
        mean_value_V50 = np.sum(V50_carrier)/measurements
        mean_value_V75 = np.sum(V75_carrier)/measurements
        mean_value_V105 = np.sum(V105_carrier)/measurements
        mean_value_V140 = np.sum(V140_carrier)/measurements
        print(mean_value_V140)
        rms_value_15 = 0
        rms_value_30 = 0
        rms_value_50 = 0
        rms_value_75 = 0
        rms_value_105 = 0
        rms_value_140 = 0
        for i in range(measurements):
            rms_value_15 += np.sqrt((V15_carrier[i]-mean_value_V15)**2)
            rms_value_30 += np.sqrt((V30_carrier[i]-mean_value_V30)**2)
            rms_value_50 += np.sqrt((V50_carrier[i]-mean_value_V50)**2)
            rms_value_75 += np.sqrt((V75_carrier[i]-mean_value_V75)**2)
            rms_value_105 += np.sqrt((V105_carrier[i]-mean_value_V105)**2)
            rms_value_140 += np.sqrt((V140_carrier[i]-mean_value_V140)**2)
        rms_value_15 = rms_value_15/measurements
        rms_value_30 = rms_value_30/measurements
        rms_value_50 = rms_value_50/measurements
        rms_value_75 = rms_value_75/measurements
        rms_value_105 = rms_value_105/measurements
        rms_value_140 = rms_value_140/measurements
        print(rms_value_140)
        

        param_bounds = ([0],[1])
        x_points = [(xedges[l+1]+xedges[l])/2.0 for l in range(len(xedges)-1)] 
     #    plt.errorbar(x_points, mean_values_15, yerr=rms_values_15, fmt='')
     #    plt.errorbar(x_points, mean_values_50, yerr=rms_values_50, fmt='')
     #    plt.errorbar(x_points, mean_values_140, yerr=rms_values_140, fmt='')
     #ax = photon_dataframe.plot.hist(bins=100, alpha=0.5)


#        print('REAL V='+str('%.2f'%V_values[j]))
#        print('RETRIEVED VISIBILITY for n='+str(n)+' photons')
#        print('BW=15nm:'+str(ufloat(popt15[0],perr15)))
#        print('BW=30nm:'+str(ufloat(popt30[0],perr30)))
#        print('BW=50nm:'+str(ufloat(popt50[0],perr50)))
#        print('BW=75nm:'+str(ufloat(popt75[0],perr75)))
#        print('BW=105nm:'+str(ufloat(popt105[0],perr105)))
#        print('BW=140nm:'+str(ufloat(popt140[0],perr140)))
#        
        Recovered_V_n_BW[j][n].extend((mean_value_V15,mean_value_V30,mean_value_V50,mean_value_V75,mean_value_V105,mean_value_V140))
        Uncertainty_V_n_BW[j][n].extend((rms_value_15[0],rms_value_30[0],rms_value_50[0],rms_value_75[0],rms_value_105[0],rms_value_140[0]))

#        plt.show()

colors = ['blue','red','green','darkorange','darkgrey']
plt.rcParams['axes.labelsize'] = 24
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['font.weight']='bold'
plt.rcParams['font.size'] = 20
plt.rcParams['legend.handlelength'] = 2 #Added by me
plt.rcParams['legend.markerscale'] = 2 #Added by me
fig2,ax2 = plt.subplots(figsize=(11,8))
markertypes=['.','+','s']
for i in range(len(V_values)):
    ax2.axhline(y=V_values[i], xmin=0, xmax=1,color= colors[i],linestyle='--', linewidth='3')
    for n in range(len(photon_numbers)): 
        plt.errorbar(bandwidth_values,Recovered_V_n_BW[i][n],yerr=Uncertainty_V_n_BW[i][n],fmt='.',marker=markertypes[n], markersize='10',color = colors[i])
    #plt.scatter(wavelengths,carrier_measured_values[i],marker='x',s=100,color = colors[i], label='Measured visibility V='+str('%.2f'%visibility_values[i]))
ax2.set_xlim(10,150)
ax2.set_ylim(0.0,1.0)
ax2.set_xlabel('Filter Bandwidth (nm)')
ax2.set_ylabel(r'Visibility $| \gamma_{12} |$')
#-----------------LEGEND------------------------
f = lambda m,c: plt.plot([],[],marker=m, color=c, ls="none")[0]
lo = plt.scatter(random(10), random(10), marker='.', color='black')
li = plt.scatter(random(10), random(10), marker='+', color='black')
ll = plt.scatter(random(10), random(10), marker='s', color='black')
plt.legend((lo,li,ll), (r'$2$x$10^4$ ph',r'$1$x$10^5$ ph',r'$1$x$10^6$ ph'))
#------------------------------------------
fig2.savefig('Raw_data_2/Retrieved_visibility_combined_new.jpg',dpi=1200)

markertypes=['.','.','.']
for n in range(len(photon_numbers)):
    fig,ax = plt.subplots(figsize=(11,8))
    for i in range(len(V_values)):
        ax.axhline(y=V_values[i], xmin=0, xmax=1,color= colors[i],linestyle='--', linewidth='3')
        plt.errorbar(bandwidth_values,Recovered_V_n_BW[i][n],yerr=Uncertainty_V_n_BW[i][n],fmt='.',marker=markertypes[n], markersize='10',color = colors[i])
        #plt.scatter(wavelengths,carrier_measured_values[i],marker='x',s=100,color = colors[i], label='Measured visibility V='+str('%.2f'%visibility_values[i]))
    ax.set_xlim(10,150)
    ax.set_ylim(0.0,1.0)
    ax.set_xlabel('Filter Bandwidth (nm)')
    ax.set_ylabel(r'Visibility $| \gamma_{12} |$')
    fig.savefig('Raw_data_2/Retrieved_visibility_'+str(photon_numbers[n])+'_new.jpg',dpi=1200)
