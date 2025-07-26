# -*- coding: utf-8 -*-
"""
Created on Thursday 29 December 2020

@author: Elif
"""
import numpy as np
import matplotlib as matplotlib
from matplotlib import pyplot as plt
import os

import fromlist_loop
import branchDirection_loop

def frequency_of_followup_check(activityType):
    frequency = activityType.sum()/len(activityType)
    return frequency


def count_distance(distance, N):
    distance_occurance = np.zeros(2*N-1) 
    for p_idx, p in enumerate(range(-(N-2),N-1,1)):
        distance_occurance[p_idx] = (distance==p).sum()/len(distance)
    return distance_occurance


def count_last_neuron_id(last_neuron_id, simulationspercombination, N):
    occurance = np.zeros(N-1) # encoded pattern length = 7
    for p in range(N-1):
        occurance[p] = (last_neuron_id==p).sum()/simulationspercombination
    return occurance


def computeEfficacity(working_folder, name_file, simulationspercombination, c, N):
    nextActivityDictionary = {x: [] for x in (np.arange(1,N))}

    eff = np.zeros(simulationspercombination)
    last_neuron_id = np.zeros(simulationspercombination)
    nextActivation = np.zeros(simulationspercombination)
    followingActivityType = np.zeros(simulationspercombination)
    distance = np.zeros(simulationspercombination)
    duration = np.zeros(simulationspercombination)
    occurance_of_distance = 0
    mean_distance = 0
    std_distance = 0

    for n in range(simulationspercombination):

        name_file_trial = name_file +"_%d"%(int(n))
        working_txt_patternActivity = os.path.join(working_folder, name_file_trial+'.txt')
        row_id, last_neuron, followingActivityType[n], nextActivation[n], distance[n] = fromlist_loop.analyseActiveUnits(working_txt_patternActivity, c)
        last_neuron_id[n] = last_neuron 
        eff[n] = last_neuron/6

        if followingActivityType[n]>0.:
            nextActivityDictionary[int(last_neuron)].append(distance[n])
            
    occurance = count_last_neuron_id(last_neuron_id, simulationspercombination,N)
    frequency_of_followup = frequency_of_followup_check(followingActivityType)
    mean_duration = np.mean(duration[np.isfinite(duration)])
    std_duration = np.std(duration[np.isfinite(duration)])

    if np.any(followingActivityType > 0.):
        occurance_of_distance = count_distance(distance[np.where(followingActivityType>0)], N)
        mean_distance = np.mean(distance[np.where(followingActivityType>0)])
        std_distance = np.std(distance[np.where(followingActivityType>0)])
    
    return np.mean(eff), occurance, frequency_of_followup, occurance_of_distance, nextActivityDictionary, mean_distance, std_distance


def last_neruon_id_matrix(working_folder, name_file, simulationspercombination, c):

    last_neuron_id = np.zeros(simulationspercombination)

    for n in range(simulationspercombination):

        name_file_trial = name_file +"_%d"%(int(n))
        working_txt_patternActivity = os.path.join(working_folder, name_file_trial+'.txt')
        _, last_neuron, _, _, _ = fromlist_loop.analyseActiveUnits(working_txt_patternActivity, c)
        last_neuron_id[n] = last_neuron 

    return last_neuron_id


def next_neruon_id_matrix(working_folder, name_file, simulationspercombination, c):

    next_neuron_id = np.zeros(simulationspercombination)
    # distance = np.zeros(simulationspercombination)

    for n in range(simulationspercombination):

        name_file_trial = name_file +"_%d"%(int(n))
        working_txt_patternActivity = os.path.join(working_folder, name_file_trial+'.txt')
        _, _, _, next_neuron_id[n], _ = fromlist_loop.analyseActiveUnits(working_txt_patternActivity, c)

    return next_neuron_id


def branch_neruon_id_matrix(working_folder, name_file, simulationspercombination, c, N, branchUnit, distance):
    
    branch_direction = np.zeros(simulationspercombination)
    branch_direction_after_deactivation = np.zeros(simulationspercombination)


    for n in range(simulationspercombination):
        name_file_trial = name_file +"_%d"%(int(n))
        working_txt_patternActivity = os.path.join(working_folder, name_file_trial+'.txt')
        branch_direction[n], branch_direction_after_deactivation[n] = branchDirection_loop.analyseBranch(working_txt_patternActivity, c, branchUnit)
    
    branch_frequency = branchDirection_loop.calculateBranchFrequency(branch_direction, simulationspercombination, N, branchUnit, distance)
    branch_frequency_after_deactivation = branchDirection_loop.calculateBranchFrequency(branch_direction_after_deactivation, simulationspercombination, N, branchUnit, distance)


    return branch_direction, branch_frequency, branch_direction_after_deactivation, branch_frequency_after_deactivation

def calculateBranchFromPattern(pattern, branchPattern, patternLegend):
    index_of_pattern = patternLegend.index(pattern)
    branch = branchPattern[index_of_pattern]
    return branch


def branch_wrt_pattern(pattern_before, pattern_after, simulationspercombination, p):
    branch_before = np.zeros(simulationspercombination)
    branch_after = np.zeros(simulationspercombination)

    for n in range(simulationspercombination):
        # print("n=%d"%(n), "pattern_before=%s"%(pattern_before[n]), "pattern_after=%s"%(pattern_after[n]))÷
        # 
        branch_before[n] = calculateBranchFromPattern(pattern_before[n], p.branchPattern, p.patternLegend)
        branch_after[n] = calculateBranchFromPattern(pattern_after[n], p.branchPattern, p.patternLegend)

    return branch_before, branch_after

def pattern_id_matrix(working_folder, name_file, p):
    
    pattern_vector = []
    pattern_vector_before_deactivation = []
    for n in range(p.simulationspercombination):
        name_file_trial = name_file +"_%d"%(int(n))
        working_txt_patternActivity = os.path.join(working_folder, name_file_trial+'.txt')
        # print(working_txt_patternActivity)
        after, before = branchDirection_loop.analysePattern(working_txt_patternActivity, p)
        # print(before, after)
        pattern_vector.append(after)
        pattern_vector_before_deactivation.append(before)
    
    pattern_frequency = branchDirection_loop.calculatePatternFrequency(pattern_vector, p.simulationspercombination, p.patternLegend)
    pattern_frequency_before_deactivation = branchDirection_loop.calculatePatternFrequency(pattern_vector_before_deactivation, p.simulationspercombination, p.patternLegend)

    return pattern_vector, pattern_frequency , pattern_vector_before_deactivation, pattern_frequency_before_deactivation

def writeFile_NextActivityDetail(nextActivityDetail, working_txt_nextPatternActivity):
    with open(working_txt_nextPatternActivity, "w+") as output:
        output.write('\n'.join(str(line) for line in nextActivityDetail))
        output.close()


def create_folder(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

##########################################
####### Parameter combinations ######
##########################################
def execute(p, eta, tau, u, mu, G):
    N = p.N

    liste_mu = p.liste_mu_change

    folder_name = "Branch=%d"%(int(p.branch_type))+"_N=%d"%(int(N))+ "_G=%d"%(int(1000*G))+ "_eta=%d"%(int(1000*eta))+\
                    "_tau=%d"%(int(tau)) + "_U=%d"%(int(1000*u)) + "_IC=%d"%(int(np.sum(p.liste_units_change[0]))) 
    
    branch_vector_before_deactivation = np.zeros((len(liste_mu),p.simulationspercombination))
    branch_vector_after_deactivation = np.zeros((len(liste_mu),p.simulationspercombination))

    pattern_vector_after_deactivation = []
    pattern_vector_before_deactivation = []

                            
    current_folder = os.path.dirname(os.path.realpath(__file__))
    working_folder = os.path.join(current_folder, folder_name)

    for iidx, mu in np.ndenumerate(p.liste_mu_change):
        name_file = folder_name + "_lmu=%d"%(int(1000*mu)) 
        pattern_after, _, pattern_before, _ = pattern_id_matrix(working_folder, name_file, p)
        pattern_vector_after_deactivation.append(pattern_after)
        pattern_vector_before_deactivation.append(pattern_before)
        branch_vector_before_deactivation[iidx], branch_vector_after_deactivation[iidx], = branch_wrt_pattern(pattern_before, pattern_after, p.simulationspercombination, p)

    return pattern_vector_before_deactivation, pattern_vector_after_deactivation, branch_vector_before_deactivation, branch_vector_after_deactivation