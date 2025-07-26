import numpy as np
import analysis_sequence
import fromlist_loop


def following_activity_check(units, unitsTime, last_neuron):
    activityType = []
    activation = [0.0]
    npa = np.asarray(unitsTime)
    
    #activation.clear()
    try:
        # print('entry try')
        activationIndex = np.where(npa[:,1]!=None)[0][0]
        # print('exit try')
    except IndexError:
        # print('entry except')
        activityType = 0
        activation[0] = last_neuron
        activationIndex = 0
    else:
        # print('entry else')
        activityType = 1
        # print('len', len(units))
        activation = units[activationIndex]
    return activityType, activation, activationIndex


def consistency_check(ref, next, coactiveUnit):
    # non consistent = -1
    # pending = 0
    # consist = 1
    if next[0] == None :
        return 0 
    else:
        if coactiveUnit == 2: #### to count 2 coavtive units
            if ref[-1] == next[0]:
                return 1
            else:
                return -1
        else: #### to count 3 coavtive units
            if ref[-1] in next:
                return 1 # last unit of the reference is in the next: regular trantions
            else:
                if (len(ref) == 2) and (len(next) == 3):
                    if ref[-1] == next[1]:
                        return 1
                elif (len(ref) == 3) and (len(next) == 2):
                    if ref[-1] == next[1]:
                        return 1 # last unit of the reference is in the next: regular trantions
                else:
                    return -1


def branch_check(data, coactiveUnit, branchUnit):
    cnt = 0
    last_seq = 0
    last_cons = 10
    
    branch_direction = 0
    branch_cons = 10
    while cnt < len(data)-1:
        ref = data[cnt]
        next = data[cnt+1]
        if ref[0] == None:
            ref = data[cnt-1]
            next = next
        # here we check the consistency
        res = consistency_check(ref, next, coactiveUnit)
        if res == -1:
            last_seq = cnt
            branch_direction = 0
            break
        elif res == 0:
            branch_cons = res
        else:
            if ref[0] == branchUnit:
                if branch_cons != 0:
                    #branch_direction = max(next) - min(next)
                    branch_direction = max(next) - branchUnit
                else:
                    try:
                        #branch_direction = max(data[cnt+2]) - min(data[cnt+2]) 
                        branch_direction = max(data[cnt+2]) - branchUnit       
                    except Exception:
                        # branch_direction = max(data[cnt+1]) - min(data[cnt+1])
                        branch_direction = max(data[cnt+1]) - branchUnit 
                break
        cnt = cnt + 1
    return branch_direction

def branch_check_after_deactivation(data, coactiveUnit, branchUnit):
    # print(data)
    cnt = 0
    last_seq = 0
    last_cons = 10
    
    branch_direction = 0
    ref = [branchUnit]
    next = []
    
    for i_index, i in enumerate(data):
        if (i[0] == None) and (i_index + 1 <len(data)):
            next.append(data[i_index+1])

    if next == []:
        branch_direction = branch_check(data, coactiveUnit, branchUnit)
        # print('branch_direction 1', branch_direction)
    else:
        print(next)
        branch_direction = next[0][-1] - ref[-1]
        # print('branch_direction 2', branch_direction)

        if branch_direction <= 0:
            branch_direction = 0
        elif branch_direction >= 4:
            branch_direction = 4
        else:
            branch_direction = 1
    
    # print('branch_direction 3', branch_direction)

    return branch_direction
    # if next == []:
        # branch_direction = branch_check(data, coactiveUnit, branchUnit)
    # else:
    #     branch_direction = next[0][-1] - branchUnit
  
    # return branch_direction÷


def pattern_check_after_deactivation(data, list_pattern_units, list_patterns):
    cnt = 0
    next = []
    pattern_units = []
    first_pattern_index = 0
    next_activation_indexes = []
    for i_index, i in enumerate(data):
        if (i[0] == None) and (i_index + 1 <len(data)):
            next.append(data[i_index+1])
            next_activation_indexes.append(i_index+1)        
    # print(data)
    # print(next)
    if next == []:
        pattern_units = [None] # this lines is to avoid the error when the pattern is not found. it is changed from data[2] to [None]
    elif len(next[0]) == 2:
        pattern_units = next[0]
    else:
        first_pattern_index = next_activation_indexes[0]+1
               
        if first_pattern_index +1 < len(data):
            pattern_units = [None]
            if (data[first_pattern_index] == [None]):
                pattern_units = [None]
            elif len(data[first_pattern_index]) == 3:
                pattern_units = data[first_pattern_index+1]
            else:
                pattern_units = data[first_pattern_index]
        elif first_pattern_index +1 == len(data):
            pattern_units = data[first_pattern_index]
        else:
            pattern_units = [None]      
    try:        
        pattern = list_patterns[list_pattern_units.index(pattern_units)]
    except:
        pattern_units = [None]
        pattern = list_patterns[list_pattern_units.index(pattern_units)]
        
    return pattern

def pattern_check_before_deactivation(activeUnits, row_id, list_pattern_units, list_patterns):
    pattern = []
    
    if len(activeUnits[row_id]) == 2:
        pattern_units = activeUnits[row_id]
        pattern = list_patterns[list_pattern_units.index(pattern_units)]
    elif len(activeUnits[row_id]) == 1:
        pattern_units = activeUnits[row_id-1]
        pattern = list_patterns[list_pattern_units.index(pattern_units)]

    return pattern

def analyseBranch(working_txt_patternActivity, coactiveUnit, branchUnit):
    activeUnits = []
    activeUnitsTimes = []
    followingActivityType = []
    nextActivation = []
    
    with open(working_txt_patternActivity, "r") as file:
        for line in file:
            f1, f2 = analysis_sequence.parse_line(line)
            activeUnits.append(f1)
            activeUnitsTimes.append(f2)
    file.close()

    branch_direction = branch_check(activeUnits, coactiveUnit, branchUnit)
    # row_id, _ = fromlist_loop.efficiency_check(activeUnits, coactiveUnit)

    branch_direction_after_deactivation = branch_check_after_deactivation(activeUnits, coactiveUnit, branchUnit)

    return branch_direction, branch_direction_after_deactivation

def analysePattern(working_txt_patternActivity,  p):
    activeUnits = []
    activeUnitsTimes = []
    followingActivityType = []
    nextActivation = []
    
    with open(working_txt_patternActivity, "r") as file:
        for line in file:
            f1, f2 = analysis_sequence.parse_line(line)
            activeUnits.append(f1)
            activeUnitsTimes.append(f2)
    file.close()
    
    row_id, _ = fromlist_loop.efficiency_check(activeUnits, 3)

    pattern_before_deactivation = pattern_check_before_deactivation(activeUnits,row_id, p.patternUnit, p.patternLegend)

    pattern_after_deactivation = pattern_check_after_deactivation(activeUnits,  p.patternUnit, p.patternLegend)
    
    return pattern_after_deactivation, pattern_before_deactivation

def calculateBranchFrequency(branch_direction_vector, simulationspercombination, N, branchUnit, distance):
    
    # backwards, stops, forward_1, forward_2
    branch_frequency = np.zeros(len(distance),) # encoded pattern length = 7

    for p in range(len(distance)):
        branch_frequency[p] = (branch_direction_vector==distance[p]).sum()/simulationspercombination
    return branch_frequency

def calculatePatternFrequency(pattern_vector, simulationspercombination, list_patterns):
    
    pattern_frequency = np.zeros(len(list_patterns),) # encoded pattern length = 7

    for p in range(len(list_patterns)):
        if simulationspercombination:
            pattern_frequency[p] = (pattern_vector.count(list_patterns[p]))/simulationspercombination
        else:
            pattern_frequency[p] = 0

    return pattern_frequency
