import numpy as np
import matplotlib.pyplot as plt

#Import the template data
K5 = np.loadtxt('MiniProjectAllData/Template/keck_k5.txt')
K5_wavelength = K5[:,0]
K5_flux = K5[:,1]
K5_sigma = K5[:,2]

#Plotting the template data
plt.figure(figsize=(1920/144, 720/144), dpi=144)
plt.plot(K5_wavelength,K5_flux)
plt.ylabel(r'$F_\lambda$')
plt.xlabel(r'$\lambda (\AA)$')
plt.grid()
plt.title('Keck K5 Template Spectrum')
plt.savefig('MiniProjectAllData\Template\keck_k5_raw.png')

#loading in the GS2000 data and plotting
for i in range(0,13):
    #Loading in the observational data (MJD and phase) for the 13 spectra 
    spectrum_data = np.loadtxt('MiniProjectAllData\spectrum_data.txt')
    #Loading each spectrum itself
    if i<9:
        GS2000 = np.loadtxt('MiniProjectAllData\GS2000\keck_gs2000_0'+str(i+1)+'.txt')
    else:
        GS2000 = np.loadtxt('MiniProjectAllData\GS2000\keck_gs2000_'+str(i+1)+'.txt')
    
    #Plotting each spectrum
    plt.figure(figsize=(1920/144, 720/144), dpi=144)
    plt.plot(GS2000[:,0],GS2000[:,1])
    plt.title('GS2000 '+str(i+1)+r', $MJD = $'+str(spectrum_data[i,1])+r', $\phi = $'+str(spectrum_data[i,2]))
    plt.xlabel(r'$\lambda (\AA)$')
    plt.ylabel(r'$F_\lambda$')    
    plt.grid()
    
    #Saving each spectrum
    if i<9:
        plt.savefig('MiniProjectAllData\GS2000\keck_gs2000_subtraction_0'+str(i+1)+'_raw.png')
    else:
        plt.savefig('MiniProjectAllData\GS2000\keck_gs2000_subtraction_'+ str(i+1)+'_raw.png')