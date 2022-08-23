import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import csv
import pandas as pd
from uncertainties import ufloat
from numpy.random import random

def sinc(x):
    return np.sin(x)/x

def sinc2(x):
    return (np.sin(x)/x)**2

def interference_pattern(position,visibility):
    a = 0.0005
    d = 0.001
    f = 0.1
    wavelength = 550e-9
    x1 = np.pi*a/(wavelength*f)*(position)
    x2 = 2*np.pi*d/(wavelength*f)*(position)
    return sinc2(x1)*(1+visibility*np.cos(x2))/(1+visibility) #Normalized to 1


def calculate_rms_size(visibility,L,uncertainty):
    a = 0.0005
    d = 0.001
    f = 0.1
    wavelength = 550e-9
    sigma = wavelength*L*np.sqrt(np.log(1/np.sqrt(visibility)))/(np.pi*d)
    delta_sigma = wavelength*L/(np.pi*d)*uncertainty/(4*visibility*np.sqrt(1/np.sqrt(visibility)))
    return sigma,delta_sigma

photon_numbers = np.array(['2e4','1e5','1e6'])
carrier_total_values= []
carrier_total_errors = []
for number in photon_numbers:
    visibility_values = np.array([0.1008,0.2303,0.4728,0.6928,0.9542])
    carrier_measured_values = [] 
    carrier_error_values = [] 
    for i in range(len(visibility_values)):
        detection_file = 'detection_photons_V_'+str('%02d'%(visibility_values[i]*100))+'_'+number+'particle.txt'
        detection_dataframe = pd.read_csv(detection_file, delim_whitespace = True)
        #detection_dataframe.columns=['z_bunch','t_bunch','Ekin_bunch']
        x_edges = np.array(detection_dataframe[detection_dataframe.columns[0]].astype(float).tolist())
        nm5 = np.array(detection_dataframe[detection_dataframe.columns[1]].astype(float).tolist())
        nm10 = np.array(detection_dataframe[detection_dataframe.columns[2]].astype(float).tolist())
        nm30 = np.array(detection_dataframe[detection_dataframe.columns[3]].astype(float).tolist())
        nm50 = np.array(detection_dataframe[detection_dataframe.columns[4]].astype(float).tolist())
        nm75 = np.array(detection_dataframe[detection_dataframe.columns[5]].astype(float).tolist())
        nm100 = np.array(detection_dataframe[detection_dataframe.columns[6]].astype(float).tolist())
        nm150 = np.array(detection_dataframe[detection_dataframe.columns[7]].astype(float).tolist())
        #nm500 = np.array(detection_dataframe[detection_dataframe.columns[8]].astype(float).tolist())
        
        #Normalize the counts to a value of 1 at the maximum
        nm5 = nm5/np.amax(nm5)
        nm10 = nm10/np.amax(nm10)
        nm30 = nm30/np.amax(nm30)
        nm50 = nm50/np.amax(nm50)
        nm75 = nm75/np.amax(nm75)
        nm100 = nm100/np.amax(nm100)
        nm150 = nm150/np.amax(nm150)
        #nm500 = nm500/np.amax(nm500)
        
        param_bounds = ([0],[1])
        x_points = x_edges + (x_edges[1]-x_edges[0])/2.0
        popt5,pcov5 = curve_fit(interference_pattern,x_points,nm5,p0=(0.7),bounds=param_bounds)
        popt10,pcov10 = curve_fit(interference_pattern,x_points,nm10,p0=(0.7),bounds=param_bounds)
        popt30,pcov30 = curve_fit(interference_pattern,x_points,nm30,p0=(0.7),bounds=param_bounds)
        popt50,pcov50 = curve_fit(interference_pattern,x_points,nm50,p0=(0.7),bounds=param_bounds)
        popt75,pcov75 = curve_fit(interference_pattern,x_points,nm75,p0=(0.7),bounds=param_bounds)
        popt100,pcov100 = curve_fit(interference_pattern,x_points,nm100,p0=(0.7),bounds=param_bounds)
        popt150,pcov150 = curve_fit(interference_pattern,x_points,nm150,p0=(0.7),bounds=param_bounds)
        #popt500,pcov500 = curve_fit(interference_pattern,x_points,nm500,p0=(0.7),bounds=param_bounds)
        
        perr5 = np.sqrt(np.diag(pcov5))[0]
        perr10 = np.sqrt(np.diag(pcov10))[0]
        perr30 = np.sqrt(np.diag(pcov30))[0]
        perr50 = np.sqrt(np.diag(pcov50))[0]
        perr75 = np.sqrt(np.diag(pcov75))[0]
        perr100 = np.sqrt(np.diag(pcov100))[0]
        perr150 = np.sqrt(np.diag(pcov150))[0]
        #perr500 = np.sqrt(np.diag(pcov500))[0]
        
        print('RETRIEVED VISIBILITY:')
        print('REAL V='+str(visibility_values[i]))
        print('BW=5nm:'+str(ufloat(popt5[0],perr5)))
        print('BW=10nm:'+str(ufloat(popt10[0],perr10)))
        print('BW=30nm:'+str(ufloat(popt30[0],perr30)))
        print('BW=50nm:'+str(ufloat(popt50[0],perr50)))
        print('BW=75nm:'+str(ufloat(popt75[0],perr75)))
        print('BW=100nm:'+str(ufloat(popt100[0],perr100)))
        print('BW=150nm:'+str(ufloat(popt150[0],perr150)))
        #print('BW=500nm:'+str(ufloat(popt500[0],perr500)))
        
        wavelengths = np.array([5.0,15.0,30.0,50.0,75.0,100.0,145.0])#,500.0])
        values = np.array([popt5[0],popt10[0],popt30[0],popt50[0],popt75[0],popt100[0],popt150[0]])#,popt500[0]])
        fit_errors = np.array([perr5,perr10,perr30,perr50,perr75,perr100,perr150])#,perr500])
        
        #plt.scatter(x_points,nm20)
        #plt.plot(x_points,interference_pattern(x_points,*popt20),color='red')
        
        fig, axs = plt.subplots(3, 2)
        fig.suptitle('V='+str('%.2f'%visibility_values[i]))
        #plt.title('V='+str('%.2f'%i))
        axs[0, 0].plot(x_points*1e3, interference_pattern(x_points,*popt5),color='red')
        axs[0, 0].scatter(x_points*1e3,nm5,s=2,zorder=5)
        axs[0, 0].set_title('BW = 5nm')
        axs[0, 1].plot(x_points*1e3, interference_pattern(x_points,*popt10),color='red',label='Fit')
        axs[0, 1].scatter(x_points*1e3,nm10,s=2,zorder=5,label='Measured')
        axs[0, 1].set_title('BW = 10nm')
        axs[1, 0].plot(x_points*1e3, interference_pattern(x_points,*popt30),color='red')
        axs[1, 0].scatter(x_points*1e3,nm30,s=2,zorder=5)
        axs[1, 0].set_title('BW = 30nm')
        axs[1, 1].plot(x_points*1e3, interference_pattern(x_points,*popt50),color='red')
        axs[1, 1].scatter(x_points*1e3,nm50,s=2,zorder=5)
        axs[1, 1].set_title('BW = 50nm')
        axs[2, 0].plot(x_points*1e3, interference_pattern(x_points,*popt75),color='red')
        axs[2, 0].scatter(x_points*1e3,nm75,s=2,zorder=5)
        axs[2, 0].set_title('BW = 75nm')
        axs[2, 1].plot(x_points*1e3, interference_pattern(x_points,*popt100),color='red')
        axs[2, 1].scatter(x_points*1e3,nm100,s=2,zorder=5)
        axs[2, 1].set_title('BW = 100nm')
        for ax in axs.flat:
            ax.set(xlabel='Position [mm]', ylabel='Intensity [a.u.]')
        # Hide x labels and tick labels for top plots and y ticks for right plots.
        for ax in axs.flat:
            ax.label_outer()
        #plt.legend()
        plt.tight_layout()
        fig.savefig('Fits_V='+str('%.2f'%visibility_values[i])+'.pdf')
        plt.close()

        carrier_measured_values.append(values)
        carrier_error_values.append(fit_errors)
        print('\n')
    carrier_total_values.append(carrier_measured_values)
    carrier_total_errors.append(carrier_error_values)



#-------------COMBINED PLOTTING------------------------
colors = ['blue','red','green','darkorange','darkgrey']
plt.rcParams['axes.labelsize'] = 24
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['font.weight']='bold'
plt.rcParams['font.size'] = 20
plt.rcParams['legend.handlelength'] = 2 #Added by me
plt.rcParams['legend.markerscale'] = 2 #Added by me


fig2,ax2 = plt.subplots(figsize=(11,8))
markertypes = ["x","+","s"]
for j in range(len(photon_numbers)):

    for i in range(len(carrier_measured_values)):
        #plt.errorbar(wavelengths,carrier_total_values[j][i],yerr=carrier_total_errors[j][i],fmt='.',marker=markertypes[j],markersize='10',color = colors[i], label='Measured visibility V='+str('%.2f'%visibility_values[i]))
        plt.scatter(wavelengths,carrier_total_values[j][i],marker=markertypes[j],s=100,color = colors[i], label='Measured visibility V='+str('%.2f'%visibility_values[i]))
#        ax2.axhline(y=visibility_values[i], xmin=0, xmax=1, color=colors[i],linestyle='--', linewidth='3')
    
for i in range(len(carrier_measured_values)):
    ax2.axhline(y=visibility_values[i], xmin=0, xmax=1, color=colors[i],linestyle='--', linewidth='3')
ax2.set_xlim(10,150)
ax2.set_ylim(0.0,1.0)
ax2.set_xlabel('Filter Bandwidth (nm)')
ax2.set_ylabel(r'Visibility $| \gamma_{12} |$')
#-----------------LEGEND------------------------
f = lambda m,c: plt.plot([],[],marker=m, color=c, ls="none")[0]
lo = plt.scatter(random(10), random(10), marker='x', color='black')
li = plt.scatter(random(10), random(10), marker='+', color='black')
ll = plt.scatter(random(10), random(10), marker='s', color='black')
plt.legend((lo,li,ll), (r'$2$x$10^4$ ph',r'$1$x$10^5$ ph',r'$1$x$10^6$ ph'))
#------------------------------------------
fig2.savefig('Retrieved_visibility_Ji_Gwang_combined.jpg', dpi=1200)
plt.show()


print('RETRIEVED RMS BEAM SIZE:')
L_values = [0.04,0.05,0.07,0.1,0.3]
carrier_size = []
carrier_error = []
for i in range(len(carrier_measured_values)):
    print('------V='+str(visibility_values[i]))
    sizes_x = np.array([])
    errors_x = np.array([])
    L = L_values[i]
    for j in range(len(carrier_measured_values[i])):
        rms_x,error_x = calculate_rms_size(carrier_measured_values[i][j],L,carrier_error_values[i][j])
        sizes_x = np.append(sizes_x,rms_x)
        errors_x = np.append(errors_x,error_x)
        print(ufloat(rms_x,error_x))
    carrier_size.append(sizes_x)
    carrier_error.append(errors_x)
    print('\n')


print('RETRIEVED AVERAGE RMS BEAM SIZE (BV<=100nm):')
for k in range(len(carrier_size)):
    average_size = 0.0
    error_average_size = 0.0
    for number in range(4): #Only for BW of 100nm or smaller
        average_size += carrier_size[k][number]
        error_average_size += carrier_error[k][number]**2
    average_size =  average_size/4.0
    error_average_size = np.sqrt(error_average_size)/4
    print('------V='+str(visibility_values[k]))
    print(ufloat(average_size,error_average_size))

