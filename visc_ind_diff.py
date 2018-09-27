# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 15:37:12 2018

Direct biofilm simulation (Kelvin-Voigt 1D chain)

@author: dyanni3
"""

#%% imports
import numpy as np
from matplotlib import pyplot as plt
from scipy import io as sio

#%% define constants and set-up
E = 1
eta = 1
T = 20000 #total num time steps
lambda_act = 1e-2 #rate of birth/death events
N = 10000 #length of biofilm
sigma_0 = 10

# positions of events
j = np.random.randint(0,N,size=T)

#nature of events (birth = 1, death = -1)
s = np.random.choice([-1,1], size=T)

#timing of events
dt = np.random.exponential(scale = 1/lambda_act, size = T)
#dt = (1/lambda_act)*np.ones(T)
t = np.cumsum(dt)

sigma = np.zeros(shape=(T,N))
epsilon = np.zeros(shape=(T,N))

#%% run the simulation
for k in range(T):
	sigma[k,:j[k]] += -sigma_0*s[k]
	sigma[k,j[k]:] += sigma_0*s[k]
	epsilon[k,:] = (1 / (1 + dt[k]*E/eta)) * (epsilon[k-1,:] + (dt[k]/eta)*sigma[k,:])
	
#%% plot results
plt.plot(t,epsilon[:,int(N/2)]**2)
plt.plot(t,np.average(epsilon**2,axis=1))

#%% run many trials
etas = [.1,1,10,100,1000]
plt.figure(figsize=(8,8))
for i,eta in enumerate(etas):
	M = 100 #num trials
	times = np.zeros((M,T))
	eps2 = np.zeros((M,T))
	for m in range(M):
		if m%10==0:
			print(m)
		j = np.random.randint(0,N,size=T)
		s = np.random.choice([-1,1], size=T)
		dt = np.random.exponential(scale = 1/lambda_act, size = T)
		t = np.cumsum(dt)
		sigma = np.zeros(shape=(T,N))
		epsilon = np.zeros(shape=(T,N))
		for k in range(T):
			sigma[k,:j[k]] = sigma[k-1,:j[k]] - sigma_0*s[k]
			sigma[k,j[k]:] = sigma[k-1,j[k]:] + sigma_0*s[k]
			epsilon[k,:] = (1 / (1 + dt[k]*E/eta)) * (epsilon[k-1,:] + (dt[k]/eta)*sigma[k,:])
		eps2[m,:] = np.average(epsilon**2,axis=1)
		times[m,:] = t
	plt.loglog(np.average(times,axis=0),np.average(eps2,axis=0),label=eta)
	sio.savemat("diff_eta%d"%i,{"eps2":np.average(eps2,axis=0)})
	sio.savemat("diff_eta_times%d"%i,{"times":np.average(times,axis=0)})
plt.legend()

#%% run many trials
etas = np.array([.1,100,1000])
lambda_acts = 1e-6*((E/etas)**1) 
plt.figure(figsize=(8,8))
plt.xlabel("Time",fontsize=16)
plt.ylabel(r"$\langle\Delta\,x^2\rangle$",fontsize=16)
for i,eta in enumerate(etas):
	M = 500 #num trials
	times = np.zeros((M,T))
	eps2 = np.zeros((M,T))
	for m in range(M):
		if m%10==0:
			print(m)
		j = np.random.randint(0,N,size=T)
		s = np.random.choice([-1,1], size=T)
		dt = np.random.exponential(scale = 1/lambda_acts[i], size = T)
		t = np.cumsum(dt)
		sigma = np.zeros(shape=(T,N))
		epsilon = np.zeros(shape=(T,N))
		for k in range(T):
			sigma[k,:j[k]] = sigma[k-1,:j[k]] - sigma_0*s[k]
			sigma[k,j[k]:] = sigma[k-1,j[k]:] + sigma_0*s[k]
			epsilon[k,:] = (1 / (1 + dt[k]*E/eta)) * (epsilon[k-1,:] + (dt[k]/eta)*sigma[k,:])
		eps2[m,:] = np.average(epsilon**2,axis=1)
		times[m,:] = t
	plt.loglog(np.average(times,axis=0),np.average(eps2,axis=0),marker='o',label=r"$\eta=$ %s"%str(eta))
	sio.savemat("diff_eta_long_passive_MSD%d"%i,{"eps2":np.average(eps2,axis=0)})
	sio.savemat("diff_eta_long_passive_times%d"%i,{"times":np.average(times,axis=0)})
plt.legend()
#%% plotting D vs \eta

plt.figure(figsize=(8,8))
plt.xlabel(r"$\eta$",fontsize=24)
plt.ylabel("D",fontsize=24)
for i in range(6):
    if(i==0 or i==5):
        plt.loglog(fake_etas,fake_Ds[i],marker='o',label=r"$\alpha=$ %s"%str(alphas[i]),lw=6,markersize=12)
    else:
        plt.loglog(fake_etas,fake_Ds[i],marker='o',label=r"$\alpha=$ %s"%str(alphas[i]))
plt.legend(fontsize=16)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#%%
etas = np.array([.1,1,10,100,1000])
for j,alpha in enumerate(alphas):
	lambda_acts = 1e-2  *((E/etas)**alpha) 
	plt.figure(figsize=(8,8))
	plt.xlabel("Time",fontsize=16)
	plt.ylabel(r"$\langle\Delta\,x^2\rangle$",fontsize=16)
	for i,eta in enumerate(etas):
		M = 1000 #num trials
		times = np.zeros((M,T))
		eps2 = np.zeros((M,T))
		for m in range(M):
			if m%10==0:
				print(m)
			j = np.random.randint(0,N,size=T)
			s = np.random.choice([-1,1], size=T)
			dt = np.random.exponential(scale = 1/lambda_acts[i], size = T)
			t = np.cumsum(dt)
			sigma = np.zeros(shape=(T,N))
			epsilon = np.zeros(shape=(T,N))
			for k in range(T):
				sigma[k,:j[k]] = sigma[k-1,:j[k]] - sigma_0*s[k]
				sigma[k,j[k]:] = sigma[k-1,j[k]:] + sigma_0*s[k]
				epsilon[k,:] = (1 / (1 + dt[k]*E/eta)) * (epsilon[k-1,:] + (dt[k]/eta)*sigma[k,:])
			eps2[m,:] = np.average(epsilon**2,axis=1)
			times[m,:] = t
		plt.loglog(np.average(times,axis=0),np.average(eps2,axis=0),marker='o',label=r"$\eta=$ %s"%str(eta))
		sio.savemat("diff_eta_MSD%s_alpha%s"%(str(eta),str(alpha)),{"eps2":np.average(eps2,axis=0)})
		sio.savemat("diff_eta_times%s_alpha%s"%(str(eta),str(alpha)),{"times":np.average(times,axis=0)})
	plt.legend()
