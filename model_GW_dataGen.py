import pylab
import csv
import numpy as np
import matplotlib.pyplot as pp
import warnings
warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")
import lal
from pycbc.waveform import get_td_waveform
from pycbc import waveform
import time
from tqdm import tqdm
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d
import glob
import os

#################
#function to find index of closest element in list
def closest(lst, K):
    return min(range(len(lst)), key = lambda i: abs(lst[i]-K))

#################
#prepare time array for saving with less resolution, s.t. time array does not loose resolution at small times
def timeArray(N_tot):
    d_time = (1-10**-1)/N_tot
    time_array = np.arange(-1, -10**-1, d_time)
    d_time = (10**-1-10**-2)/N_tot
    time_array = np.append(time_array,np.arange(-10**-1, -10**-2, d_time))
    d_time = (10**-2-10**-3)/N_tot
    time_array = np.append(time_array,np.arange(-10**-2, -10**-3, d_time))
    d_time = (10**-3-10**-4)/N_tot
    time_array = np.append(time_array,np.arange(-10**-3, -10**-4, d_time))
    time_array = np.append(time_array,-10**-4)
    return time_array

#################
#function to read and store rho and cs matrices
def readMatrix(namefile):
    list_to_fill = []
    with open(namefile, mode='r') as file:
        reader = csv.reader(file, delimiter = '\t')
        for row in reader:
            list_to_fill.append(row)
    return list_to_fill

#################
#read mass and tidal def array
def readMassTidal(eos,column_mass,column_tidal):
    #rescaling factor to be in agreement with PyCBC units
    KK = (2.997 * 10.0**5)**2/(6.67 * 10.0**-20 * 1.988 * 10.0**30)
    mass_array = []
    tidal_array = []
    with open(eos, mode='r') as file: 
        reader = csv.reader(file, delimiter =' ') 
     
        for row in reader: 
            mass_array.append(float(row[column_mass]))
            tidal_array.append(float(row[column_tidal])*KK**5/
                                float(row[column_mass])**5)
    return mass_array, tidal_array

#################
#task in main loop (integration for each EOS)
def task_EOS_loop(eos):
    
    #get specifics of EOS from the file name
    parts = eos.split('_')
    #info on low energy eos
    if parts[1] == 'ap4':
        low_eos = 0
    if parts[1] == 'sly':
        low_eos = 1
    #info on cosmological constant value
    Lambda_value = float(parts[2])
    #info on QCD modification (correspond to a specific line in c_s and rho file
    z_value = int(parts[len(parts)-1].split('.')[0])
      
    #reset ranges of grid to standard
    M_min = 0.9 
    M_max = 2.5
    deltaM = 0.1
    MAX = int((M_max-M_min)/deltaM)
    
    #Interpolate tidal def vs M
    column_index_mass = 1
    column_index_tidal = 3
    mass_values, tidal_values = readMassTidal(eos,column_index_mass,column_index_tidal)
 
    index_max = mass_values.index(max(mass_values))
    #if M_max is smaller than 2.5 I restrict the grid
    if max(mass_values) < M_max:
        M_max = max(mass_values)
        MAX = int((M_max-M_min)/deltaM)
        
    interp_tidal_mass = interp1d(mass_values[0:index_max+1],tidal_values[0:index_max+1])
    
    list_to_save = [] #list to save waveform for the whole grid
    #loop on m1
    for i in tqdm(range(0,MAX+1), desc='Loop on m1 for eos={}'.format(int(z_value)), leave=True):
        #loop on m2
        for j in range(i,MAX+1):
            ###
            m1 = round(0.9 + i * deltaM,1)      
            m2 = round(0.9 + j * deltaM,1)
            l1 = float(interp_tidal_mass(m1))
            l2 = float(interp_tidal_mass(m2))
            freq = 80 #By a brief analysis this is a good choice to have signal ~1s
        
            approximant = "TEOBResumS"

            hp, hc = get_td_waveform(approximant = approximant,
                                     mass1 = m1, 
                                     mass2 = m2,
                                     lambda1 = l1, 
                                     lambda2 =l2,
                                     #smaller time step increases precision and avoid errors
                                     delta_t = 1.0 / 4096, 
                                     f_lower = freq)

            hp, hc = hp.trim_zeros(), hc.trim_zeros()

            #amp = waveform.utils.amplitude_from_polarizations(hp, hc)
            #f = waveform.utils.frequency_from_polarizations(hp, hc)

            #CHECK for t=-1
            #index = closest(hp.sample_times, -1)
            #print(hp.sample_times[0],hp.sample_times[index])
            #print(index, hp.sample_times[index],f[1],f[index]) #CHECK for best model
        
            #Extrapolate only 1002 points on last second of signal:
            time_array = timeArray(250)
            interp_waveform = CubicSpline(hp.sample_times,hp)
            new_interped_waveform = interp_waveform(time_array)
            
            list_to_save_temp = np.insert(new_interped_waveform, 0, [z_value, low_eos, Lambda_value,
                                                                m1, m2, l1, l2], axis=0)
            list_to_save.append(list_to_save_temp)
            
    return list_to_save