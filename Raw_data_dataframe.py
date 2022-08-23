import pandas as pd

bandwidth_values = [14,29,50,75,105,140]
V_values = [0.100855,0.230338,0.472799,0.692774,0.960037]
for j in V_values:
    detected_photons = pd.DataFrame({})
    carrier_raw_data = []
    for i in bandwidth_values:
        lst = []
        with open('Raw_data_3/'+str('%d'%(i))+'BW_'+str('%2f'%j)+'V.txt', 'r') as infile:
            for line in infile:
                line = line.strip('\n')
                lst.append(float(line))
        
        detected_photons[str('%d'%i)] = lst
        #carrier_raw_data.append(lst)
        print(1)
    print('hallo!')
#    detected_photons = pd.DataFrame(
#            {str('%d'%(bandwidth_values[0]*1e9))+'nm': carrier_raw_data[0],
#        str('%d'%(bandwidth_values[1]*1e9))+'nm': carrier_raw_data[1],
#        str('%d'%(bandwidth_values[2]*1e9))+'nm': carrier_raw_data[2],
#        str('%d'%(bandwidth_values[3]*1e9))+'nm': carrier_raw_data[3],
#        str('%d'%(bandwidth_values[4]*1e9))+'nm': carrier_raw_data[4],
#        str('%d'%(bandwidth_values[5]*1e9))+'nm': carrier_raw_data[5],
#        })
    newfile = 'Raw_data_3/Photons_V_'+str('%02d'%(j*100))+'_raw_data.txt'   #Write the projection in x into a txt file, this will simulate an experimental readout of a screen with 2um pixels and 1e5 incident photons
    detected_photons.to_csv(newfile, sep='\t', float_format='%.7f',index=False)
