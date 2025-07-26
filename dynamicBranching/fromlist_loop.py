import numpy as np
import analysis_sequence


def following_activity_check(units, unitsTime, last_neuron ):
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

def sequence_check(data, coactiveUnit):
    cnt = 0
    last_seq = 0
    last_cons = 10
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
            break
        elif res == 0:
            last_cons = res
        else:
            last_seq = cnt+1
        cnt = cnt + 1
    return last_seq

def efficiency_check(data, coactiveUnit):
    row_id = sequence_check(data, coactiveUnit)
    last_neuron = data[row_id][-1]
    if last_neuron == None:
        row_id = row_id -1
        last_neuron = data[row_id][-1]
        return row_id, last_neuron
    else:
        return row_id, last_neuron

def seqence_duration_check(last_time_row, last_neuron):
    a = [a for a in last_time_row if a != None]
    sequence_duration = a[0] / last_neuron
    # for aidx, a in enumerate(last_time_row):
    #     if a != None:
    #         if aidx == 0:
    #             sequence_duration = a /last_neuron
    #         else:
    #             sequence_duration = np.inf
    return sequence_duration

def check_distance_nextActivation(nextActivation, last_neuron):
    #print(nextActivation, type(nextActivation))
    #if nextActivation != last_neuron:
    nextActivatedUnit = [0]
    if len(nextActivation) == 1:
            distance = nextActivation[0] - last_neuron
            nextActivatedUnit[0] = nextActivation[0]
    else:
        for x in nextActivation:
            if x not in [last_neuron]:
                distance = x - last_neuron
                nextActivatedUnit[0] = x
#    else:
#        distance = nextActivation[0] - last_neuron
    return distance, nextActivatedUnit

def analyseActiveUnits(working_txt_patternActivity, coactiveUnit):
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

    row_id, last_neuron = efficiency_check(activeUnits, coactiveUnit)

    # print('inital', last_neuron)
    #sequence_duration = seqence_duration_check(activeUnitsTimes[row_id], last_neuron)
    ########################################################################################################
###################### AFTER THE INITIAL CHECK I CONTINUE WITH THE FOLLOWING ACTIVITY ##########################
    ########################################################################################################

    followingActivityType, nextActivation, nextActivationIndex = following_activity_check(activeUnits[row_id+1:], activeUnitsTimes[row_id+1:], last_neuron)
    followingActivationRowId = nextActivationIndex + row_id+1

####### continue without rethreating the distance
    # if followingActivityType == 0:
    #     nextActivatedUnit = nextActivation[0]
    #     distance = 0
    # else:
    #     if len(nextActivation) == 1:
    #         distance = nextActivation[0] - last_neuron
    #         nextActivatedUnit = nextActivation[0]
    #     else: 
    #         for x in nextActivation:
    #             if x not in [last_neuron]:
    #                 distance = x - last_neuron
    #                 nextActivatedUnit = x
            
    if followingActivityType == 0:
        nextActivatedUnit = nextActivation[0]
        distance = 0
    else:
        while followingActivityType == 1:
            if len(nextActivation) == 1:
                distance = nextActivation[0] - last_neuron
                nextActivatedUnit = nextActivation[0]
            else: 
                for x in nextActivation:
                    if x not in [last_neuron]:
                        distance = x - last_neuron
                        nextActivatedUnit = x
            if distance == 0: 
                # HERE I COUNT FOR LOOPING. IF len(nextActivation) == 1: distance = nextActivation[0] - last_neuron GIVES SELF ACTIVATION
                # I RE-CHECK THE EFFICIENCY STARTING FROM THE LAST ACTIVATED UNIT. 
                row_id, last_neuron = efficiency_check(activeUnits[followingActivationRowId:], coactiveUnit)
                row_id = followingActivationRowId+row_id
                followingActivityType, nextActivation, nextActivationIndex = following_activity_check(activeUnits[row_id+1:], activeUnitsTimes[row_id+1:], last_neuron)
                followingActivationRowId = nextActivationIndex + row_id+1
                #print('here:', row_id, last_neuron)
            else:
                break
            
    
    return row_id, last_neuron, followingActivityType, nextActivatedUnit, distance 