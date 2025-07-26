import os
import numpy as np

class parameters_linear:
    def __init__(self):
        self.simulationspercombination = 5 # simulations per combination
        
        self.branch_type = 1 # for Case 1: linear

        self.distance = np.array([0, 0., 0.])
        self.N=10
        self.branch_unit = [3] # the unit from which the sequence bifurcates 

        self.liste_eta = np.array([0.04])
        self.liste_rho = np.array([1.2])
        self.liste_tau = np.array([300]) 
        self.liste_G = np.array([0.600])
        self.liste_mu = np.array([0.4])
        self.liste_mu_change = np.array([0.4])
        self.liste_units_change = [[4, 5]] # the unit from which the sequence bifurcates
        self.I=0.0

        self.simTime = 2000
        self.dt = 0.1 

        self.patternLegend = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'None']
        self.unitLegend = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']

        self.patternUnit = [[0,1], [1,2], [2,3], [3,4], [4,5], [5,6], [6,7], [7, 8], [8, 9], [None]]
        self.branchPattern = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]   

class parameters_branching:
    def __init__(self):
        self.simulationspercombination = 5 # simulations per combination
        
        # self.branch_type = 1 # for Case 1: linear
        self.branch_type = 2 # for Case 2: branching
        # self.branch_type = 3 # for Case 3: cycling

        self.distance = np.array([0, 1., 4.])
        self.N=10
        self.branch_unit = [3] # the unit from which the sequence bifurcates 

        self.liste_eta = np.array([0.04])
        self.liste_rho = np.array([1.2])
        self.liste_tau = np.array([300]) 
        self.liste_G = np.array([0.600])
        self.liste_mu = np.array([0.4])
        self.liste_mu_change = np.array([0.5])
        self.liste_units_change = [[4, 5]] # the unit from which the sequence bifurcates
        self.I=0.0

        self.simTime = 2000
        self.dt = 0.1 

        self.patternLegend = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'None']
        self.unitLegend = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']

        self.patternUnit = [[0,1], [1,2], [2,3], [3,4], [4,5], [5,6], [3,7], [7, 8], [8, 9], [None]]
        self.branchPattern = [0, 0, 0, 1, 1, 1, 4, 4, 4, 0]   


class parameters_cycling:
    def __init__(self):
        self.simulationspercombination = 5 # simulations per combination
        
        self.branch_type = 3 # for Case 3: cycling

        self.distance = np.array([0, 1., 4.])
        self.N=10
        self.branch_unit = [3] # the unit from which the sequence bifurcates 

        self.liste_eta = np.array([0.04])
        self.liste_rho = np.array([1.2])
        self.liste_tau = np.array([300]) 
        self.liste_G = np.array([0.600])
        self.liste_mu = np.array([0.4])
        self.liste_mu_change = np.array([0.5])
        self.liste_units_change = [[4, 5]] # gain of these units will change by liste_mu_change
        self.I=0.0

        self.simTime = 2000
        self.dt = 0.1 

        self.patternLegend = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'None']
        self.unitLegend = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']

        self.patternUnit = [[0,1], [1,2], [2,3], [3,4], [4,5], [5,6], [3,7], [7, 8], [8, 9], [6, 9], [None]]
        self.branchPattern = [0, 0, 0, 1, 1, 1, 4, 4, 4, 5, 0]