# -*- coding: utf-8 -*-
"""
@author: Elif Koksal Ersoz
The system of N units is simulated using the Euler-Maruyama method for dt =0.01.
Units above the threshold 0.5 are considered to be active. To faciitate the analysis, solutions are filtered (smoothened).
Passages from the threshold are recorded for furhter analysis.

"""
import numpy as np
import os
from scipy.signal import savgol_filter
from matplotlib import colors as mcolors

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

################################################
##### Synaptic plasticity in E-E matrices ######
################################################

def create_folder(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def writeActivityPattern(working_txt_patternActivity, activityPattern,activityTime): 
    with open(working_txt_patternActivity, "w") as output:
        indexActivityTime = 0
        for p in activityPattern[1:]:
            output.write('{!s};{!s}'.format(p, activityTime[indexActivityTime]))
            output.write("\n")
            indexActivityTime += 1
        output.close()

def patternChoose(branchType):
    if branchType == 1: 
        pattern = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0],
                            [1, 1, 0, 0, 0, 0, 0, 0, 0],
                            [0, 1, 1, 0, 0, 0, 0, 0, 0],
                            [0, 0, 1, 1, 0, 0, 0, 0, 0],
                            [0, 0, 0, 1, 1, 0, 0, 0, 0],
                            [0, 0, 0, 0, 1, 1, 0, 0, 0],
                            [0, 0, 0, 0, 0, 1, 1, 0, 0],
                            [0, 0, 0, 0, 0, 0, 1, 1, 0],
                            [0, 0, 0, 0, 0, 0, 0, 1, 1],
                            [0, 0, 0, 0, 0, 0, 0, 0, 1]])
    if branchType == 2: 
        pattern = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0],
                            [1, 1, 0, 0, 0, 0, 0, 0, 0],
                            [0, 1, 1, 0, 0, 0, 0, 0, 0],
                            [0, 0, 1, 1, 0, 0, 1, 0, 0],
                            [0, 0, 0, 1, 1, 0, 0, 0, 0],
                            [0, 0, 0, 0, 1, 1, 0, 0, 0],
                            [0, 0, 0, 0, 0, 1, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 1, 1, 0],
                            [0, 0, 0, 0, 0, 0, 0, 1, 1],
                            [0, 0, 0, 0, 0, 0, 0, 0, 1]])
    elif branchType == 3:
        pattern = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 1, 1, 0, 0, 1, 0, 0, 0],
                            [0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 1, 0, 0, 0, 1],
                            [0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0, 1, 1]])
    return(pattern)

def connection(N, pattern):
	wCons = np.zeros((N, N))

	for h in range(N):
		for j in range(N):
			for p in range(np.size(pattern, axis=1)):
				wCons[h, j] = wCons[h, j]+(pattern[h, p])*(pattern[j, p])
				
	return wCons

def connectionInhibition(branchType, N, G):
    Jinh = G*np.ones((N,N))
    Jinh[3,3] = branchType*G 
    return Jinh


def find_in_iterable(x, iterable):
    for i, item in enumerate(iterable):
        if item == x:
            return i
    return None

def model(t, x, I, U, tau, mu, J, Jinh):
    # inhibition by the excitatory population 
    # Excitatiory neuron
    N = J.shape[0]
    dxdt = np.zeros((1, 2*N))

    S = x[N:2*N]*x[0:N]
    i=1
    
    for i in range(N):
	    dxdt[0][i] = x[i]*(1-x[i])*(-mu[i]*x[i] -I - np.dot(Jinh[i,0:N], x[0:N])+ np.dot(J[i,0:N], S.flatten("F")))
        
    dxdt[0][N:2*N] = (1 - x[N:2*N])/tau - U*x[N:2*N]*x[0:N]

    dxdt = dxdt.reshape((1, 2*N))
    
    return dxdt

def dW(delta_t, N):
    return np.sqrt(delta_t)*(1-2*np.random.random_sample((1,N)))

def simulation(p, eta, Tau, u, mu, G):
	
    # asymétrique
    pattern = patternChoose(p.branch_type)
    # Jmax-Hebb rule
    Jmax = connection(p.N, pattern)

    # Jinh modified wrt to the pattern
    Jinh = connectionInhibition(p.branch_type, p.N, G)

    t0 = 0.0
    ts  = np.arange(t0, p.simTime, p.dt)
    t = 1.0
    
	# Initial conditions
    x0 = np.empty(2*p.N)
    x0[0:p.N] = pattern[:,0].conj().transpose() # E neurons
    x0[p.N:2*p.N] = np.ones((1, p.N)) # S values
	
    xComplete = np.zeros((int(ts.size),2*p.N))	
    xComplete[0] = x0.copy()	

    #####################
    # Simulate the system
    # ##################		
    Winc = np.empty(2*p.N)

    for simCount in range(1, ts.size):
        
        t = (simCount - 1)*p.dt

        Winc[0:p.N]  = dW(p.dt, p.N)
        Winc[p.N:2*p.N] = np.zeros(p.N)
        
        xComplete[simCount] = xComplete[simCount-1] + model(t,xComplete[simCount-1], p.I, u, Tau, mu, Jmax, Jinh) * p.dt + eta * Winc
        
        # xi values are bounded by [0,1]
        for i in range(p.N):
            if xComplete[simCount][i]<0:
                xComplete[simCount][i]=-xComplete[simCount][i]
            elif xComplete[simCount][i]>1:
                xComplete[simCount][i]=2-xComplete[simCount][i]
    
    return ts, xComplete

def save_activity_pattern(ts, xComplete, p, eta, Tau, u, mu, G, n):
    
    # define folders, file names	
    current_folder = os.path.dirname(os.path.realpath(__file__))
    # create a folder for the current simulation
    folder_name = "Branch=%d"%(int(p.branch_type))+"_N=%d"%(int(p.N))+ "_G=%d"%(int(1000*G))+ "_eta=%d"%(int(1000*eta))+\
            "_tau=%d"%(int(Tau)) + "_U=%d"%(int(1000*u))+ "_IC=%d"%(int(np.sum(p.liste_units_change[0])))  
    working_folder = os.path.join(current_folder, folder_name)
    create_folder(working_folder)
        
    name_file = folder_name + "_lmu=%d"%(int(1000*mu[p.branch_unit[0]+1])) +"_%d"%(int(n))
    
    working_txt_patternActivity= os.path.join(working_folder, name_file+'.txt')
    # List to record neuron activities
    liste_neurones = [[] for _ in range(p.N)]
    activityPattern = [[0,1]]
    # pattern durations
    activityTime = []
    # find active neurons which are above 0.5 after applying a filter
    xCompleteFilterred = np.zeros((int(ts.size),p.N))
    for i in range(p.N):
        xCompleteFilterred[:,i] = savgol_filter(xComplete[:,i], 1001, 2)
    
    for simCount in range(ts.size):
        liste_neurones.append([])
        for  countActive, x_i in np.ndenumerate(xCompleteFilterred[simCount]):
            if x_i > 0.5:
                liste_neurones[-1].extend(list(countActive))

        if liste_neurones[-1] != liste_neurones[-2]:
            if len(liste_neurones[-1]) < len(liste_neurones[-2]): 
                activityTime.append([simCount*p.dt,''])
                activityPattern.append(liste_neurones[-1])
            elif len(liste_neurones[-1]) > len(liste_neurones[-2]): 
                activityTime.append(['',simCount*p.dt])
                activityPattern.append(liste_neurones[-1])
    
    writeActivityPattern(working_txt_patternActivity, activityPattern, activityTime)

