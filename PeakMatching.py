import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import scipy.constants as const

###Calculates and performs doppler shift on a template spectrum###
###Velocity is given in c###
def dopplershift(template, velocity):
    doppler_shift = np.sqrt((1+velocity*1000/const.c)/(1-velocity*1000/const.c))*template[:,0]
    template_fn = UnivariateSpline(doppler_shift,template[:,1],ext=1,s=0)
    return template_fn  #returns a function

###Calculates the scaling factor given a set of measurements
###and their fitted function
def scalingfactor(measurements,function):
    numerator = np.sum(measurements[:,1]*function(measurements[:,0])/measurements[:,2]**2)
    denominator = np.sum((function(measurements[:,0])**2)/(measurements[:,2]**2))
    if denominator == 0:
        denominator = np.nan
    return numerator/denominator
    
###Return chi-squared values given a set of measurements 
###and their fitted function
def chisquared(measurements,function):
    x_i = measurements[:,0]
    y_i = measurements[:,1]
    sigma_i = measurements[:,2]
    A = scalingfactor(measurements,function)
    P = function(x_i)
    return np.sum(((y_i-A*P)/(sigma_i))**2)

template = np.loadtxt('MiniProjectAllData\Template\SubtractedTemplate\keck_k5_subtracted.txt')  
spectrum_data = np.loadtxt('MiniProjectAllData\spectrum_data.txt')  
v = np.linspace(-0.005*const.c/1000,0.005*const.c/1000,1000)
chi_squared = np.zeros(len(v))
velocity = np.zeros((13,2))

for i in range(1,14):
    if i < 10:
        continuum = np.loadtxt('MiniProjectAllData/GS2000/ContinuumSubtraction/keck_gs2000_subtraction_0'+str(i)+'.txt')
    if i >= 10:
        continuum = np.loadtxt('MiniProjectAllData/GS2000/ContinuumSubtraction/keck_gs2000_subtraction_'+str(i)+'.txt')

    for j in range(0, len(v)):
        template_fn = dopplershift(template,v[j])
        chi_squared[j] = chisquared(continuum,template_fn)
    print(i)
    
    shifted = chi_squared-(np.min(chi_squared)+1)
    shifted_chi_squared_fn = UnivariateSpline(v,shifted,ext=1,s=0)
    root = shifted_chi_squared_fn.roots()
    sigma_v = root[1]-root[0]
    print(sigma_v)
    
    velocity[i-1,0] = v[np.where(chi_squared == np.min(chi_squared))]
    velocity[i-1,1] = sigma_v
    
    plt.figure(figsize=(1920/144, 1080/144), dpi=144)
    plt.plot(v,chi_squared)
    plt.xlabel(r'$v(kms^{-1})$')
    plt.ylabel(r'$\chi^2$')
    plt.axvline(velocity[i-1,0], color='r', linestyle='dashed', linewidth=1,label="v = "+str(round(velocity[i-1,0], 1))+r'$\pm$'+str(round(sigma_v,1))+r'$kms^{-1}$')
    plt.title('GS2000 '+str(i)+r', $MJD = $'+str(spectrum_data[i-1,1])+r', $\phi = $'+str(spectrum_data[i-1,2]))
    plt.grid()
    plt.legend()
    if i < 10:
        plt.savefig('MiniProjectAllData/PeakMatching/GS2000_ChiSquared_0'+str(i)+'.png')  
    if i >= 10:
        plt.savefig('MiniProjectAllData/PeakMatching/GS2000_ChiSquared_'+str(i)+'.png')  
        
    f,axarr = plt.subplots(2,1,figsize=(1920/144, 1080/144), dpi=144)
    matched_fn = dopplershift(template,velocity[i-1,0])
    axarr[0].plot(continuum[:,0],matched_fn(continuum[:,0]))
    axarr[1].plot(continuum[:,0],continuum[:,1])
    if i < 10:
        plt.savefig('MiniProjectAllData/PeakMatching/Fitted/GS2000_ChiSquared_0'+str(i)+'.png')  
    if i >= 10:
        plt.savefig('MiniProjectAllData/PeakMatching/Fitted/GS2000_ChiSquared_'+str(i)+'.png')  
    
    
np.savetxt('MiniProjectAllData/PeakMatching/velocities.txt',velocity)