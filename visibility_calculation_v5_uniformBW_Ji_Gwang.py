import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random
import csv
from scipy.stats import gaussian_kde
from matplotlib import gridspec

def sinc(x):
    return np.sin(x)/x

def sinc2(x):
    return (np.sin(x)/x)**2


def visibility_monochromatic(wavelength,sigma,d,L):
    return np.exp(-2*(np.pi*sigma*d/wavelength/L)**2)

def individual_interference(position,wavelength,initial_pos,a,d,L,f):
    x1 = np.pi*a/(wavelength*f)*(position+f*initial_pos/L)
#    x1 = np.pi*a/(wavelength*f)*(position)
    x2 = 2*np.pi*d/(wavelength*f)*(position+f*initial_pos/L)
    return sinc2(x1)*(1+np.cos(x2))/2.0 #Normalized to 1
 

def initial_montecarlo(central,sigma,limits):
    for i in range(100000):
        candidate = central + ((limits[1]-limits[0])*random.random() + limits[0])
        treshold = random.random()
        if treshold <= np.exp(-(candidate-central)**2/(2*sigma**2)):
            return candidate
        else:
            continue

def initial_uniform(central,sigma,limits):
    return central + ((limits[1]-limits[0])*random.random() + limits[0])


n= 1001
low_lim = -0.0002
high_lim = 0.0002
positions = np.linspace(low_lim,high_lim,n) #m, positions in final_screen

##BANDWIDTH (uniform shaped distribution)
#bandwidth_rms = 100.0e-9 #20.0e-9 #m
bandwidth_center = 550e-9 #m, green light
#bandwidth_limited = [-50e-9,50e-9] #m, limits for monte-carlo


#BUNCH (gaussian shaped distribution)
bunch_rms_x = 7.5e-6 #2e-6 #m, bunch rms size
bunch_center_x = 0.0 #m, centered at 0
bunch_limited_x = [-80e-6,80e-6] #m, limits for monte-carlo
bunch_rms_y = 2.3e-6 #2e-6 #m, bunch rms size
bunch_center_y = 0.0 #m, centered at 0
bunch_limited_y = [-20e-6,20e-6] #m, limits for monte-carlo

#GEOMETRY
L = 0.05 #m, distance from source to double slit
f = 0.1 #m, distance from double slit to screen
a = 0.0005 #m, width of slits
d = 0.001 #m, slit_separation

expected_V = visibility_monochromatic(bandwidth_center,bunch_rms_x,d,L)
print('For monochromatic light we would expect: '+str(expected_V))

bandwidth_values = [5.0e-9,10.0e-9,30.0e-9,50.0e-9,75.0e-9,100e-9,150e-9]#,200.0e-9,300e-9,400e-9,500e-9]
carrier_projection_x = []
for i in bandwidth_values:
    #BANDWIDTH (uniform shaped distribution)
    bandwidth_rms = i #20.0e-9 #m
    bandwidth_center = 550e-9 #m, green light
    bandwidth_limited = [-i/2,i/2] #m, limits for monte-carlo

    n_particles = 20000
    wavelengths = np.array([])
    initial_positions_x = np.array([])
    initial_positions_y = np.array([])
    final_pos_x = np.array([])
    accumulator = np.zeros(n)
    counter= 0
    for particle in range(n_particles):
        initial_pos_x = initial_montecarlo(bunch_center_x,bunch_rms_x,bunch_limited_x) #We give to each photon an initial position x
        initial_pos_y = initial_montecarlo(bunch_center_y,bunch_rms_y,bunch_limited_y) #We give to each photon an initial position y
        wavelength = initial_uniform(bandwidth_center,bandwidth_rms,bandwidth_limited) #We give to each photon an initial wavelength
        #wavelength = bandwidth_center
        wavelengths = np.append(wavelengths,wavelength)
        initial_positions_x = np.append(initial_positions_x,initial_pos_x)
        initial_positions_y = np.append(initial_positions_y,initial_pos_y)
        probab_x = individual_interference(positions,wavelength,initial_pos_x,a,d,L,f)
    #    accumulator = accumulator + probab
        for j in range(100000):
            random_index = int(n*random.random())
            random_value = random.random()
            if random_value <= probab_x[random_index]:
                final_pos_x = np.append(final_pos_x,positions[random_index])
                break
            else:
                continue
        counter+=1
        #print(counter)
    
    #accumulator = accumulator/n_particles
##    xy = np.vstack([final_pos_x,initial_positions_y])
##    z = gaussian_kde(xy)(xy)
    
    
    counts, xedges,yedges = np.histogram2d(final_pos_x,initial_positions_y,bins=[200,100],range=[[low_lim, high_lim], [-0.00001, 0.00001]]) #3um pixel size asumed in x (range -0.00075m to 0.00075m divided in 300 pixels)
    projection_x = np.sum(counts,axis=1)
    projection_y = np.sum(counts,axis=0)
    carrier_projection_x.append(projection_x)    
  
    print('Histo done')
    
    fig0,ax01 = plt.subplots()
    ax01.set_xlabel(r'x [mm]')
    ax01.set_ylabel('y [mm]')
    ax01.set_title(r'Measured Intensity')
    #ax01.hist(final_pos*1e3,bins=1001,density=True)#,label='Final photon positions')
    #ax01.plot(positions*1e3,accumulator,label= r'BW = '+str('%d' % (bandwidth_rms*1e9))+r'nm and $\sigma_{x}=$'+str('%d' % (bunch_rms*1e6))+r'$\mu$m')
    #ax01.scatter(final_pos_x,initial_positions_y,c=z,s=1)
    ax01.hist2d(final_pos_x*1e3,initial_positions_y*1e3,[300,50])
    ax01.ticklabel_format(style = 'sci')
    plt.grid()
    fig0.tight_layout()
    #plt.legend()
    fig0.savefig('Bandwidth_pure_montecarlo_image_'+str('%d'%(i*1e9))+'nmBW_2e4particle_Ji_Gwang.png')
    plt.close()
  
    print('Plot1 done')

    
    fig2=plt.figure(figsize=(6, 6))
    xlim = (xedges[:-1].min(), xedges[:-1].max())
    ylim = (yedges[:-1].min(), yedges[:-1].max())
    gs = gridspec.GridSpec(2, 2, width_ratios=[3,1], height_ratios=[3,1])
    ax = plt.subplot(gs[0,0])
    ax.set_title(r'Measured Intensity')
    axr = plt.subplot(gs[0,1], sharey=ax)
    axb = plt.subplot(gs[1,0], sharex=ax)
    ax.hist2d(final_pos_x*1e3,initial_positions_y*1e3,[300,100])
    axr.plot(projection_y, yedges[:-1]*1e3)
    axb.plot(xedges[:-1]*1e3, projection_x)
    #axr.set_ylim(ylim)
    #axb.set_xlim(xlim)
    axb.set_xlabel('x [mm]')
    axr.set_ylabel('y [mm]')
    ax.tick_params(
        axis='both',        # changes apply to the x-axis
        which='both',       # both major and minor ticks are affected
        bottom=False,       # ticks along the bottom edge are off
        left=False,         # ticks along the top edge are off
        labelbottom=False,  # labels along the bottom edge are off
        labelleft=False)    # labels along the bottom edge are off
    axb.tick_params(
        axis='y',        # changes apply to the x-axis
        which='both',       # both major and minor ticks are affected
        left=False,         # ticks along the top edge are off
        labelleft=False)    # labels along the bottom edge are off
    axr.tick_params(
        axis='x',        # changes apply to the x-axis
        which='both',       # both major and minor ticks are affected
        bottom=False,         # ticks along the top edge are off
        labelbottom=False)    # labels along the bottom edge are off
    axr.yaxis.set_label_position("right")
    axr.yaxis.tick_right()
    fig2.tight_layout()
    fig2.savefig('Bandwidth_pure_montecarlo_'+str('%d' % (bandwidth_rms*1e9))+'nm_'+str('%d' % (bunch_rms_x*1e6))+'um_xrms'+str('%d' % (bunch_rms_y*1e6))+'um_yrms_projections_2e4particle_Ji_Gwang.pdf')
    plt.close()
    
    print('Plot2 done')
    
    fig3,ax31 = plt.subplots()
    ax31.set_xlabel(r'$\lambda$ [nm]')
    ax31.set_ylabel('Photon # [normalized density]')
    ax31.set_title(r'Wavelength distribution')
    #ax31.hist(bins[:-1], bins, weights=counts)#,label='Wavelength Distribution')
    ax31.hist(wavelengths*1e9,bins=501,density=True, label='Wavelength Distribution')
    ax31.set_xlim((bandwidth_center+bandwidth_limited[0])*1e9,(bandwidth_center+bandwidth_limited[1])*1e9)
    ax31.ticklabel_format(style = 'sci')
    plt.grid()
    fig3.tight_layout()
    plt.legend()
    fig3.savefig('Wavelength_distribution_'+str('%d' % (bandwidth_rms*1e9))+'nmBW_2e4particle_Ji_Gwang.pdf')
    plt.close()

    print('Plot3 done')


detected_photons = pd.DataFrame(
        {'Pos. [m]':xedges[:-1],
     str('%d'%(bandwidth_values[0]*1e9))+'nm': carrier_projection_x[0],
     str('%d'%(bandwidth_values[1]*1e9))+'nm': carrier_projection_x[1],
     str('%d'%(bandwidth_values[2]*1e9))+'nm': carrier_projection_x[2],
     str('%d'%(bandwidth_values[3]*1e9))+'nm': carrier_projection_x[3],
     str('%d'%(bandwidth_values[4]*1e9))+'nm': carrier_projection_x[4],
     str('%d'%(bandwidth_values[5]*1e9))+'nm': carrier_projection_x[5],
     str('%d'%(bandwidth_values[6]*1e9))+'nm': carrier_projection_x[6],
    })
newfile = 'detection_photons_V_'+str('%02d'%(expected_V*100))+'_2e4particle.txt'   #Write the projection in x into a txt file, this will simulate an experimental readout of a screen with 2um pixels and 1e5 incident photons
detected_photons.to_csv(newfile, sep='\t', float_format='%.6f',index=False)
#with open(newfile, 'w') as f:
#    writer = csv.writer(f, delimiter='\t')
#    for element in range(len(projection_x)):
#        writer.writerow(['%.6f' % xedges[element], '%d' % carrier_projection_x[0][element]])
#    f.close()



fig4,ax41 = plt.subplots()
ax41.set_xlabel(r'$x_0$ [um]')
ax41.set_ylabel('Photon # [normalized density]')
ax41.set_title(r'Initial position distribution')
#ax31.hist(bins[:-1], bins, weights=counts)#,label='Wavelength Distribution')
ax41.hist(initial_positions_x*1e6,bins=501,density=True, label='Position Distribution')
ax41.set_xlim((bunch_center_x+bunch_limited_x[0])*1e6,(bunch_center_x+bunch_limited_x[1])*1e6)
ax41.ticklabel_format(style = 'sci')
plt.grid()
fig4.tight_layout()
plt.legend()
fig4.savefig('Initial_position_distribution_'+str('%d' % (bunch_rms_x*1e6))+'um_2e4particle.pdf')
plt.show()
