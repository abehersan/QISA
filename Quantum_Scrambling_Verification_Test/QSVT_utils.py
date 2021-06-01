################################################################################
# GENERAL STRUCTURING FUNCS:
################################################################################

# func to get qubit inums of upper and lower half of the circ 
def fst_n_sec_half_nums(qnum):
    fst_half =  [x for x in range(int(qnum/2))]
    sec_half = [x for x in range(int(qnum/2), qnum-1)]    
    return fst_half, sec_half
    
# find initial bell pair inums for n qubit circ, args = number of qubits (int)
def get_bell_pairs(qnum):
    initalsys = [x for x in range(qnum)][1:-1] # put first and last qubits aside for initial entanglement of the system
    pairs = [[item,initalsys[int(len(initalsys)/2):][number]] for (number,item) in enumerate(reversed(initalsys[:int(len(initalsys)/2)]))]
    pairs.append([qnum-2, qnum-1])  # append last pair with bob qubit and next neighbor
    return pairs

# func to get 2-qubit-gate pairs for scrambling operation applied to n qubit circuit, args = number of qubits (int)
def get_unitary_pairs(qnum):  
    upper_qs = [x for x in range(qnum-1)][:int(qnum/2)]
    lower_qs = [x for x in range(qnum-1)][int(qnum/2):]  
    upper_pairs = [(max(upper_qs), min(upper_qs))] + sorted([(upper_qs[i], upper_qs[i+1]) for i in range(len(upper_qs)-1)])
    lower_pairs = sorted([(lower_qs[i],lower_qs[i+1]) for i in range(len(lower_qs)-1)])                         
    lower_pairs = [(min(lower_qs),max(lower_qs))] + [x for x in reversed(lower_pairs)] 
    return upper_pairs, lower_pairs

# func to get dict with all lists plus reversed lists for qubit inum pairings, one-q-gates and two-q-gates
def get_all_gatesets_and_pairing_lists(n, one_q_gates, two_q_gates):
    '''args: 
        n (int), 
        one_q_gates, two_q_gates (lists storing gate objects)
    
    returns: 
    fst_2gate_pairs, rev_fst_2gate_pairs, sec_2gate_pairs, rev_sec_2gate_pairs, one_q_gates, rev_one_q_gates, two_q_gates, rev_two_q_gates
    '''
    qnum = 2*n +1 # num of qubits in circ
    
    # get pairings of qubit inmus for both halfs of circ
    fst_2gate_pairs, sec_2gate_pairs = get_unitary_pairs(qnum)
    
    # lists in reversed order for all items but the first one
    rev_fst_2gate_pairs = [fst_2gate_pairs[0]] + [x for x in reversed(fst_2gate_pairs[1:])]
    rev_sec_2gate_pairs = [sec_2gate_pairs[0]] + [x for x in reversed(sec_2gate_pairs[1:])]
    rev_two_q_gates =  [two_q_gates[0]] + [x for x in reversed(two_q_gates[1:])]
    rev_one_q_gates = [one_q_gates[0]] + [x for x in reversed(one_q_gates[1:])]

    # repeat gates if selection in list is smaller than num of pairs in unitary
    if len(two_q_gates) < len(fst_2gate_pairs):
        two_q_gates = (1+int(len(fst_2gate_pairs)/len(two_q_gates)))*two_q_gates
    if len(one_q_gates) < n*2:
        one_q_gates = (1+int(n*2/len(one_q_gates)))*one_q_gates 
    if len(rev_two_q_gates) < len(rev_fst_2gate_pairs):
        rev_two_q_gates = (1+int(len(rev_fst_2gate_pairs)/len(rev_two_q_gates)))*rev_two_q_gates
    if len(rev_one_q_gates) < n*2:
        rev_one_q_gates = (1+int(n*2/len(rev_one_q_gates)))*rev_one_q_gates
    
    # cut gate lists to len of pairing lists / num of qubits if they are longer
    two_q_gates = [two_q_gates[i] for i in range(len(fst_2gate_pairs))]
    one_q_gates = [one_q_gates[i] for i in range(n*2)]
    rev_two_q_gates = [rev_two_q_gates[i] for i in range(len(rev_fst_2gate_pairs))]
    rev_one_q_gates = [rev_one_q_gates[i] for i in range(n*2)]
    
    # storing all in list dict lsd
    lsd = {"upper_pairs": fst_2gate_pairs,"upper_pairs_rev": rev_fst_2gate_pairs, 
         "lower_pairs":sec_2gate_pairs, "lower_pairs_rev": rev_sec_2gate_pairs, 
         "one_q_gates":one_q_gates,"one_q_gates_rev":rev_one_q_gates, 
         "two_q_gates":two_q_gates, "two_q_gates_rev":rev_two_q_gates}
    
    return lsd

################################################################################
# GENERAL CIRCUIT OPERATION FUNCS:
################################################################################

# funcs for bell pair entanglement, args: circuit, nested list of qubit inums
def entangle_bell_pairs(circ, bell_pairs):  
    for pair in bell_pairs:
        circ.h(pair[0])
        circ.cx(pair[0],pair[1])
        
# funcs for single bell pair disentanglement
def disentangle_bell_pair(circ, pair):
    circ.cx(pair[0], pair[1])
    circ.h(pair[0])
    
# func to apply classical bob gates after base measurement of bell pair   
def apply_bob_gates(circ,ibob,measpair):
    measpair = sorted(measpair)
    circ.cz(measpair[0],ibob)
    circ.cx(measpair[1],ibob)
    
################################################################################
# FOR SCRAMBLING VERIFICATION TESTS: GATESET METHOD
################################################################################

# func to apply method "gateset" to circ, args: circ and initial lists with gates

def apply_gate_set_test(circ, two_q_gates, one_q_gates):
    
    # get n (unitary size) from qubit num
    n = int((len(circ.qubits)-1)/2) 
    
    # get list dict
    lsd = get_all_gatesets_and_pairing_lists(n, one_q_gates, two_q_gates)
    
    # and all the lists
    fst_2gate_pairs, rev_fst_2gate_pairs = lsd.get("upper_pairs"), lsd.get("upper_pairs_rev")
    sec_2gate_pairs, rev_sec_2gate_pairs = lsd.get("lower_pairs"), lsd.get("lower_pairs_rev")
    one_q_gates, rev_one_q_gates = lsd.get("one_q_gates") , lsd.get("one_q_gates_rev")
    two_q_gates, rev_two_q_gates = lsd.get("two_q_gates"), lsd.get("two_q_gates_rev")

    # apply two qubit gates for pairs in fist half of circ
    for i in range(len(fst_2gate_pairs)):
        circ.append(two_q_gates[i],fst_2gate_pairs[i])

    # apply two qubit gates for pairs in second half of circ
    for i in range(len(sec_2gate_pairs)):
        circ.append(two_q_gates[i],sec_2gate_pairs[i])

    # apply single qubit gates for all but bob qubit
    for i in range(2*n):
        circ.append(one_q_gates[i],[i])

    # apply two qubit gates for reversed pairs in fist half of circ
    for i in range(len(rev_fst_2gate_pairs)):
        circ.append(rev_two_q_gates[i],rev_fst_2gate_pairs[i])

    # apply two qubit gates for reversed pairs in second half of circ
    for i in range(len(rev_sec_2gate_pairs)):
        circ.append(rev_two_q_gates[i],rev_sec_2gate_pairs[i])
    
################################################################################
# FOR SCRAMBLING VERIFICATION TESTS: OPERATOR METHOD
################################################################################

def apply_operators(circ, operatorlist):
            
    # get n (unitary size) from qubit num
    qnum = len(circ.qubits)
    
    # inums of upper and lower half qubits
    fst_half, sec_half = fst_n_sec_half_nums(qnum)
    
    for op in operatorlist:
        circ.append(op, fst_half)
        circ.append(op.transpose(), sec_half)
     
    
################################################################################
# TUNABLE SCRAMBLING TOOLS: GENERAL
################################################################################

# func to return theta angle from parameter alpha
def get_theta_from_alpha(alpha):
    from numpy import pi
    return (alpha*pi)/2

    
################################################################################
# TUNABLE SCRAMBLING TOOLS: PULSE SCHEDULES
################################################################################

from qiskit import pulse, schedule
import numpy as np
from numpy import pi

# get control channel and cross resonance pulse implemented in CX basis gate
# for two qubit index nums and backend
from qiskit.pulse import *
from qiskit.pulse.library import *

def get_CRs(q1,q2,backend):
    inst_sched_map = backend.defaults().instruction_schedule_map 
    cx = inst_sched_map.get('cx', qubits=[q1,q2])        
    idx = 0
    while (type(cx.instructions[idx][1].channels[0]) is not ControlChannel) or \
        (type(cx.instructions[idx][1]) is not Play):
        idx += 1
    return (cx.instructions[idx][1].channels[0], cx.instructions[idx][1].pulse)

################################################################################

# The functions "get_direct_rx_sched" and "get_CR_theta_sched" are based on: 
# 'Optimized Quantum Compilation for Near-Term Algorithms with OpenPulse' arXiv:2004.11205
# Source code: https://github.com/singular-value/optimizations_via_openpulse
from pulse_compiler_helper_fns import *

# func to get pulse for single qubit rotation RX(theta)
def get_direct_rx_sched(theta, qubit_index, backend):
    '''
    args: 
            theta (angle, float),
            qubit_index (qubit index numbers, int),
            backend
            
    returns: pulse schedule        
    '''
    circ_inst_map = backend.defaults().instruction_schedule_map 
    
    x_instructions = circ_inst_map.get('x', qubits=[qubit_index]).instructions
    
    assert len(x_instructions) == 1
    
    x_samples = x_instructions[0][1].pulse.get_waveform().samples
    
    area_under_curve = sum(map(np.real, x_samples))
    
    if theta > np.pi:
        theta -= 2 * np.pi
    
    direct_rx_samples = rescale_samples(x_samples, (theta / np.pi))
    direct_rx_samplepulse = Waveform(direct_rx_samples).samples
    direct_rx_command = Play(Waveform(direct_rx_samplepulse),(backend.configuration().drive(qubit_index)))
    
    return Schedule([0, direct_rx_command])

# func to get pulse for two qubit cross resonance pulse CR(theta)
def get_CR_theta_sched(theta,control,target,backend):
    '''
    args: 
            theta (angle, float),
            control,target (qubit index numbers, int),
            backend
            
    returns: pulse schedule        
    '''
    inst_sched_map = backend.defaults().instruction_schedule_map 
    cx_instructions = inst_sched_map.get('cx', qubits=[control, target]).instructions
    xc_instructions = inst_sched_map.get('cx', qubits=[target, control]).instructions

    cr_uchan, cr_pulse = get_CRs(control,target,backend)

    inst_sched_map = backend.defaults().instruction_schedule_map 
    cx_instructions = inst_sched_map.get('cx', qubits=[control, target]).instructions
    xc_instructions = inst_sched_map.get('cx', qubits=[target, control]).instructions

    cr_control_inst = [(y.pulse,y.channel) for (x,y) in cx_instructions if type(y.channel)==ControlChannel and type(y)==Play]

    cr_drive_inst = [(y.pulse,y.channel) for (x,y) in cx_instructions if type(y.channel)==DriveChannel and type(y)==Play and type(y.pulse)==GaussianSquare]

    cr_control_inst = cr_control_inst[0]  # driving of control qubit at target's frequency
    cr_drive_inst = cr_drive_inst[0] # active cancellation tone

    flip = False

    if theta < 0:
        flip = True
        theta = -1 * theta

    if theta > 2 * np.pi:
        theta -= 2 * np.pi

    cr_uchan, cr_pulse = get_CRs(control,target,backend)

    full_area_under_curve = sum(map(np.real, cr_pulse.get_waveform().samples))
    target_area_under_curve = full_area_under_curve * (theta / (np.pi / 2))

    # CR pulse samples have gaussian rise, flattop, and then gaussian fall.
    # we want to find the start and end indices of the flattop
    flat_start = 0
    while cr_pulse.get_waveform().samples[flat_start] != cr_pulse.get_waveform().samples[flat_start + 1]:
        flat_start += 1
    assert cr_pulse.get_waveform().samples[flat_start] == cr_pulse.get_waveform().samples[flat_start + 1]

    flat_end = flat_start + 1
    while cr_pulse.get_waveform().samples[flat_end] == cr_pulse.get_waveform().samples[flat_end + 1]:
        flat_end += 1
    assert cr_pulse.get_waveform().samples[flat_end] == cr_pulse.get_waveform().samples[flat_end - 1]

    area_under_curve = sum(map(np.real, cr_pulse.get_waveform().samples[:flat_start]))
    area_under_curve += sum(map(np.real, cr_pulse.get_waveform().samples[flat_end+1:]))
    flat_duration = (target_area_under_curve - area_under_curve) / np.real(cr_pulse.get_waveform().samples[flat_start])
    flat_duration = max(0, int(flat_duration + 0.5))
    duration = len(cr_pulse.get_waveform().samples[:flat_start]) + flat_duration + len(cr_pulse.get_waveform().samples[flat_end+1:])
    if duration % 16 <= 8 and flat_duration > 8:
        flat_duration -= duration % 16
    else:
        flat_duration += 16 - (duration % 16)

    cr_drive_samples = np.concatenate([
        cr_pulse.get_waveform().samples[:flat_start],
        [cr_pulse.get_waveform().samples[flat_start]] * flat_duration,
        cr_pulse.get_waveform().samples[flat_end+1:]
    ])

    cr_control_samples = np.concatenate([
        cr_pulse.get_waveform().samples[:flat_start],
        [cr_pulse.get_waveform().samples[flat_start]] * flat_duration,
        cr_pulse.get_waveform().samples[flat_end+1:]
    ])

    assert len(cr_drive_samples) % 16 == 0
    assert len(cr_control_samples) % 16 == 0

    current_area_under_curve = sum(map(np.real, cr_control_samples))
    scaling_factor = target_area_under_curve / current_area_under_curve

    cr_drive_samples *= scaling_factor
    cr_control_samples *= scaling_factor

    cr_p_schedule = Play(Waveform(cr_drive_samples), cr_drive_inst[1]) | Play(Waveform(cr_control_samples),cr_control_inst[1])
    cr_m_schedule = Play(Waveform(-1*cr_drive_samples),cr_drive_inst[1]) | Play(Waveform(-1*cr_control_samples),cr_control_inst[1])

    if flip:
        schedule = cr_m_schedule
        schedule |= inst_sched_map.get('x', qubits=[control]) << schedule.duration
        schedule |= cr_p_schedule << schedule.duration
    else:
        schedule = cr_p_schedule
        schedule |= inst_sched_map.get('x', qubits=[control]) << schedule.duration
        schedule |= cr_m_schedule << schedule.duration
        
    return schedule

################################################################################

# get dict with control channels and cr pulses for all coupled pairs   
def get_CR_dict(theta, backend):
    crdict = dict(zip([tuple(pair) for pair in backend.configuration().coupling_map],
                      [get_CR_theta_sched(theta,pair[0],pair[1],backend) for pair in backend.configuration().coupling_map]))
    return crdict
    
################################################################################
################################################################################

    
    
    