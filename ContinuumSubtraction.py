import numpy as np
import matplotlib.pyplot as plt

h=51
###Centrally moving mean###
def CMA(f,h):
    N=len(f)
    fCMA = np.zeros(N)
    for i in range(0,int(h/2)):
        fCMA[i] = np.mean(f[:i+int(h/2)+1])
        fCMA[-i-1] = np.mean(f[-i-int(h/2)-1:])
    for i in range(int(h/2),N-(int(h/2))):
        fCMA[i] = np.mean(f[i-int(h/2):i+int(h/2)+1])
    return fCMA

###Centrally moving median###
def CMM(f,h):
    N=len(f)
    fCMM = np.zeros(N)
    for i in range(0,int(h/2)):
        fCMM[i] = np.median(f[:i+int(h/2)+1])
        fCMM[-i-1] = np.median(f[-i-int(h/2)-1:])
    for i in range(int(h/2),N-(int(h/2))):
        fCMM[i] = np.median(f[i-int(h/2):i+int(h/2)+1])
    return fCMM


###Loading in K5 template data###
K5 = np.loadtxt('MiniProjectAllData\Template\keck_k5.txt')
K5_wavelength = K5[:,0]
K5_flux = K5[:,1]
K5_sigma = K5[:,2]

###Our new template###
new_template = np.copy(K5)

#Calculating CMM on the Flux
K5_flux_CMM = CMM(K5_flux,h)

###Inputting our continuum-subtracted flux into out new template###
new_template[:,1] = K5_flux - K5_flux_CMM

np.savetxt('MiniProjectAllData\Template\SubtractedTemplate\keck_k5_subtracted.txt', new_template)

###Plotting K5###
f,axarr = plt.subplots(2,sharex = True,figsize=(1920/144, 1080/144), dpi=144)
axarr[0].plot(K5_wavelength,K5_flux)
axarr[0].plot(K5_wavelength,K5_flux_CMM, 'r-')
axarr[0].set_title('K5 Template Spectrum')
axarr[0].set_ylabel(r'$F_\lambda$')
axarr[1].plot(new_template[:,0],new_template[:,1])
axarr[1].set_xlabel(r'$\lambda (\AA)$')
plt.savefig('MiniProjectAllData\Template\keck_k5.png')

spectrum_data = np.loadtxt('MiniProjectAllData\spectrum_data.txt')
for i in range(0,13):
    if i<9:
        GS2000 = np.loadtxt('MiniProjectAllData\GS2000\keck_gs2000_0'+str(i+1)+'.txt')
    else:
        GS2000 = np.loadtxt('MiniProjectAllData\GS2000\keck_gs2000_'+str(i+1)+'.txt')
    
    GS2000_CMM = np.copy(GS2000)
    GS2000_sub = np.copy(GS2000)
    GS2000_CMM[:,1] = CMM(GS2000[:,1],h)
    GS2000_sub_flux = GS2000[:,1]-GS2000_CMM[:,1]
    GS2000_sub[:,1] = GS2000_sub_flux
    
    f,axarr = plt.subplots(2,sharex = True,figsize=(1920/144, 1080/144), dpi=144)
    axarr[0].plot(GS2000[:,0],GS2000[:,1])
    axarr[0].plot(GS2000[:,0],GS2000_CMM[:,1], 'r-' )
    axarr[0].grid(True)
    axarr[0].set_title('GS2000 '+str(i+1)+r', $MJD = $'+str(spectrum_data[i,1])+r', $\phi = $'+str(spectrum_data[i,2]))
    axarr[1].set_xlabel(r'$\lambda (\AA)$')
    axarr[0].set_ylabel(r'$F_\lambda$')    
    axarr[1].plot(GS2000[:,0],GS2000_sub[:,1])
    axarr[1].grid(True)
    
    if i<9:
        plt.savefig('MiniProjectAllData\GS2000\ContinuumSubtraction\keck_gs2000_subtraction_0'+str(i+1)+'.png')
        np.savetxt('MiniProjectAllData\GS2000\ContinuumSubtraction\keck_gs2000_subtraction_0'+str(i+1)+'.txt',GS2000_sub)
    else:
        plt.savefig('MiniProjectAllData\GS2000\ContinuumSubtraction\keck_gs2000_subtraction_'+str(i+1)+'.png')
        np.savetxt('MiniProjectAllData\GS2000\ContinuumSubtraction\keck_gs2000_subtraction_'+str(i+1)+'.txt',GS2000_sub)

    
    