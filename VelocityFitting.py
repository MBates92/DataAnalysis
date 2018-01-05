import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
import scipy.constants as const

velocity = np.loadtxt('MiniProjectAllData/PeakMatching/velocities.txt')
phase = np.loadtxt('MiniProjectAllData/spectrum_data.txt')

P = 0.3440915*24*3600
dates = phase[:,1]
sigma_P = np.mean(np.diff(dates))*24*3600

def sin_vel(phi,gamma,K_x,K_y):
    V = gamma + K_x*np.sin(2*np.pi*phi)+K_y*np.cos(2*np.pi*phi)
    #V_fn = UnivariateSpline(phi, V, ext=1)
    return V

def rad_vel_sem_amp(K_x,K_y):
    return np.sqrt(K_x**2+K_y**2)

def mass_function(q,M_c,f):
    return q**3 + 2*q**2+q-(M_c/f)

plt.figure(figsize=(1920/144, 1080/144), dpi=144)
plt.errorbar(phase[:,2],velocity[:,0],yerr = velocity[:,1],capsize = 2,fmt = 'bo',markersize=4, label = 'Data Points')
plt.xlabel(r'$\phi(rad)$')
plt.ylabel(r'$v(kms^{-1})$')
plt.grid()


init_vals = [1, 0, 1]
best_vals, covar = curve_fit(sin_vel, phase[:,2],velocity[:,0], p0=init_vals)
print(best_vals)
variance = [covar[0,0],covar[1,1],covar[2,2]]
sigma_vals = np.sqrt(variance)
print(sigma_vals)

phasespace = np.linspace(np.min(phase[:,2]),np.max(phase[:,2]),1000)
vel_fit = sin_vel(phasespace,best_vals[0],best_vals[1],best_vals[2])
plt.plot(phasespace,vel_fit,label = 'Fitted Curve')
plt.legend()
plt.title('Velocity Curve')
plt.savefig('MiniProjectAllData/PeakMatching/VelocityCurve.png')

K = rad_vel_sem_amp(best_vals[1],best_vals[2])
sigma_K = np.sqrt(((best_vals[1]*sigma_vals[1])**2+(best_vals[2]*sigma_vals[2])**2)/(best_vals[1]**2+best_vals[2]**2))

f = (P*((K*1000)**3))/(2*np.pi*const.G*2e30)
sigma_f = (3*f*sigma_K*1000)/(K*1000)

M_c = 0.7
i=90

function_space = np.linspace(-2.5,2.5,2000)
function = mass_function(function_space,M_c,f)
function_spline = UnivariateSpline(function_space,function,ext=1,s=0)
q = function_spline.roots()

M_x = M_c/q

dfdq = abs(-1*M_c*np.sin(i)**3*(3*q+1)/(q**2*(q+1)**3))
sigma_q = 1/dfdq*sigma_f
sigma_Mx = abs(M_x*sigma_q/q)

data = np.array([[best_vals[0],sigma_vals[0]],
                 [best_vals[1],sigma_vals[1]],
                 [best_vals[2],sigma_vals[2]],
                 [K,sigma_K],
                 [f,sigma_f],
                 [M_x,sigma_Mx],
                 [P,sigma_P]])
    
np.savetxt('MiniProjectAllData/data.txt',data)