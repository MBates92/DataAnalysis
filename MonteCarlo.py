import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

data = np.loadtxt('MiniProjectAllData/data.txt')

K_mean = data[3,0]
K_std = data[3,1]

P_mean = data[6,0]
P_std = data[6,1]

Mc_mean = 0.7
Mc_std = (1.0-0.4)/5.0

def gaussian(x,mean,sigma):
    p = np.exp(-((x-mean)**2)/(2*sigma**2))/(np.sqrt(2*np.pi*sigma**2))
    return p

def MCMC(N,burn_in,mean,std,function):
    theta_current = np.random.random()*mean
    std_prop = np.random.random()*std
    p_current = function(theta_current,mean,std)
    s=[]
    
    for i in range(0,N):
        theta_prop = st.norm(theta_current, std_prop).rvs()
        p_prop = gaussian(theta_prop,mean,std)
        if p_prop>=p_current:
            p_current = p_prop
            theta_current = theta_prop
        else:
            u_rnd = np.random.random()
            p_move = p_prop/p_current
            p_move = min(p_move,1)
            if u_rnd < p_move:
                p_current = p_prop
                theta_current = theta_prop
        if i >= burn_in:
            s.append(theta_current)
        print(i)
    
    s=np.array(s)
    return s

N=1000000
burn = 50000

K_samples = MCMC(N,burn,K_mean,K_std,gaussian)
plt.figure()
plt.hist(K_samples, bins = 100,normed = True)
plt.xlabel(r'$K(ms^{-1})$')
plt.ylabel(r'$N$')
plt.grid()

P_samples = MCMC(N,burn,P_mean,P_std,gaussian)
plt.figure()
plt.hist(P_samples, bins = 100,normed = True)
plt.xlabel(r'$P(s)$')
plt.ylabel(r'$N$')
plt.grid()

Mc_samples = MCMC(N,burn,Mc_mean,Mc_std,gaussian)
plt.figure()
plt.hist(Mc_samples, bins = 100,normed = True)
plt.xlabel(r'$M_c(M_odot)$')
plt.ylabel(r'$N$')
plt.grid()

N=N-burn

i_samples = np.zeros(N)
for i in range(0,N):
    i_samples[i] = np.random.random()*70 + 20
plt.figure()
plt.hist(i_samples, bins = 100, normed = True)
plt.xlabel(r'$i(\degree)$')
plt.ylabel(r'$N$')
plt.grid()

'''Mc_samples = np.zeros(N)
for i in range(0,N):
    Mc_samples[i] = np.random.random()*0.6 + 0.4
plt.figure()
plt.hist(Mc_samples, bins = 100, normed = True)
plt.xlabel(r'$M_c(M_odot)$')
plt.ylabel(r'$N$')
plt.grid()'''

f,axarr = plt.subplots(2,2,figsize=(1920/144, 1080/144), dpi=144)

axarr[0,0].hist(K_samples, bins = 100,normed = True)
axarr[0,0].set_xlabel(r'$K(kms^{-1})$')
axarr[0,0].set_ylabel(r'$N$')
axarr[0,0].grid()

axarr[0,1].hist(P_samples, bins = 100,normed = True)
axarr[0,1].set_xlabel(r'$P(s)$')
axarr[0,1].set_ylabel(r'$N$')
axarr[0,1].grid()

axarr[1,0].hist(i_samples, bins = 100, normed = True)
axarr[1,0].set_xlabel(r'$i(\degree)$')
axarr[1,0].set_ylabel(r'$N$')
axarr[1,0].grid()

axarr[1,1].hist(Mc_samples, bins = 100, normed = True)
axarr[1,1].set_xlabel(r'$M_c(M_\odot)$')
axarr[1,1].set_ylabel(r'$N$')
axarr[1,1].grid()

plt.savefig('MiniProjectAllData/prior.png')

samples = np.zeros((N,4))
samples[:,0] = K_samples
samples[:,1] = P_samples
samples[:,2] = i_samples
samples[:,3] = Mc_samples

np.savetxt('MiniProjectAllData/samples.txt',samples)