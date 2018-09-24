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
T = 10000 #total num time steps
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
	M = 1000 #num trials
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