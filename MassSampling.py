import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

samples = np.loadtxt('MiniProjectAllData/samples.txt')

K = samples[:,0]
P = samples[:,1]
inc = samples[:,2]
Mc = samples[:,3]
    
f = (P*(K*1000)**3)/(2*np.pi*const.G*2e30)

inc_rad = inc*np.pi/180.
denom = np.sin(inc_rad)**3
Mx = f/denom

print(np.min(Mx))

Mx[np.where(Mx>100)] = np.nan

mean = np.nanmean(Mx)
median = np.nanmedian(Mx)
std = np.nanstd(Mx)

x=[]

f,axarr = plt.subplots(2,1,figsize=(1920/144, 1080/144), dpi=144)

axarr[0].hist(Mx[~np.isnan(Mx)], normed = True, bins = 100)
axarr[0].plot(x,color = 'w',label = r"$\sigma = $"+str(round(std, 1))+
     r'$M_\odot$')
axarr[0].axvline(mean, color='r', linestyle='dashed', linewidth=1,
     label="Mean = "+str(round(mean, 1))+r'$M_\odot$')
axarr[0].axvline(median, color='g', linestyle='dashed', linewidth=1,
     label="Median = "+str(round(median, 1))+r'$M_\odot$')
axarr[0].axvline(3, color='b', linestyle='dashed', linewidth=1,
     label=r"$M_x = 3.0M_\odot$")
axarr[0].set_xlabel(r'$M_x(M_\odot)$')
axarr[0].set_ylabel(r'$p(M_x|K,P,i,M_c)$')
axarr[0].legend(shadow=True, fancybox=True)
axarr[0].grid()
axarr[1].hist(Mx[~np.isnan(Mx)], normed = True, bins = 100, cumulative = True)
axarr[1].axvline(mean, color='r', linestyle='dashed', linewidth=1)
axarr[1].axvline(median, color='g', linestyle='dashed', linewidth=1)
axarr[1].axvline(3, color='b', linestyle='dashed', linewidth=1)
axarr[1].set_xlabel(r'$M_x(M_\odot)$')
axarr[1].set_ylabel(r'$P(M_x|K,P,i,M_c)$')
axarr[1].grid()

plt.savefig('MiniProjectAllData/CompactMassDist.png')