import matplotlib, qiskit, pprint
from matplotlib import pyplot as plt
import numpy as np
from numpy import pi
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, Aer, IBMQ, transpile, schedule, assemble
import qiskit

# funcs for bell pair entanglement
def entangle_bell_pairs(circ, bell_pairs):
    # Define many-body-system bell pairs
    for pair in bell_pairs:
        circ.h(pair[0])
        circ.cx(pair[0],pair[1])
        
# funcs for bell pair disentanglement
def disentangle_bell_pair(circ, pair):
    circ.cx(pair[0], pair[1])
    circ.h(pair[0])
    
# func to insert bell pair measurement
def insert_bell_measurement(circ, pair):
    circ.barrier()
    circ.cx(pair[0],pair[1])
    circ.h(pair[0])
    circ.measure(pair[0],pair[0])
    circ.measure(pair[1],pair[1])
    circ.barrier()
    
# find initial bell pairs for n qubit circ
def get_bell_pairs(circ):
    n = len(circ.qubits)
    # make pairings:
    if n%2 != 0:
        initalsys = [x for x in range(n)][1:-1] # put first and last qubits aside for initial entanglement of the system
        pairs = [[item,initalsys[int(len(initalsys)/2):][number]]
                         for (number,item) in enumerate(reversed(initalsys[:int(len(initalsys)/2)]))]
    else:
        print("Find out how to deal with even number circ, Hannah!")
    # append pair bob and q before bob:   
    pairs.append([n-2, n-1]) 
    return pairs


# func to get structure of n qubit circ
def get_unitary_pairs(circ):  
    n = len(circ.qubits)
    if n%2 != 0:
        upper_qs = [x for x in range(n-1)][:int(n/2)]
        lower_qs = [x for x in range(n-1)][int(n/2):]
    else:
        print("Find out how to deal with even number circ, Hannah!")
    
    upper_pairs = sorted([(upper_qs[i],upper_qs[i+1]) 
                   for i in range(len(upper_qs)-1)] + [(upper_qs[i],upper_qs[i+2]) 
                                                       for i in range(len(upper_qs)-2)])

    lower_pairs = sorted([(lower_qs[i],lower_qs[i+1]) 
                   for i in range(len(lower_qs)-1)] + [(lower_qs[i],lower_qs[i+2]) 
                                                       for i in range(len(lower_qs)-2)])
    
    lower_pairs = [x for x in reversed(lower_pairs)]  
    return upper_pairs, lower_pairs

def apply_grover(circ,pair):
    from numpy import pi
    circ.rz(pi,pair[0])
    circ.rx(pi,pair[0])
    circ.rx(pi,pair[1])
    circ.swap(pair[0],pair[1])
    circ.rz(pi,pair[0])
    
########################################################################
# TUNABLE SCRAMBLING CIRCS
########################################################################


# func to return thet for Z rotation from Unitary parameter alpha

def get_theta_from_alpha(alpha):
    return (alpha*pi)/2

# func to apply tunable scrambling Unitary (upper half of circuit)

def tunable_srambling_U(circ, alpha):
    
    upper_pairs, _ =  get_unitary_pairs(circ)   
    theta = get_theta_from_alpha(alpha)
    qubit_inums = [x for x in range(int(len(circ.qubits)/2))]
    
    # U and U*
    for pair in upper_pairs:
        circ.rxx(-pi/4, pair[0],pair[1]) # rxx -pi/4 all
    circ.barrier(qubit_inums)
    for i in qubit_inums:
        circ.rz(-theta,i)                # rz -theta all
    for pair in upper_pairs:
        circ.rxx(-pi/4, pair[0],pair[1]) # rxx -pi/4 all
    circ.barrier(qubit_inums)
    for i in qubit_inums:
        circ.rz(theta,i)                 # rz +theta all
        
# func to apply tunable scrambling Unitary transpose (lower half of circuit)

def tunable_srambling_U_transpose(circ, alpha):
    
    _, lower_pairs =  get_unitary_pairs(circ)   
    theta = get_theta_from_alpha(alpha)
    qubit_inums = [x for x in range(int(len(circ.qubits)/2),len(circ.qubits))]
    
    # U and U*
    for pair in lower_pairs:
        circ.rxx(-pi/4, pair[0],pair[1]) # rxx -pi/4 all
    circ.barrier(qubit_inums)
    for i in qubit_inums:
        circ.rz(-theta,i)                # rz -theta all
    for pair in lower_pairs:
        circ.rxx(-pi/4, pair[0],pair[1]) # rxx -pi/4 all
    circ.barrier(qubit_inums)
    for i in qubit_inums:
        circ.rz(theta,i)                 # rz +theta all
        

########################################################################
# RESULT ANALYSIS TOOLS
########################################################################

# funcs to get distribution of numbers of exited states in measurement result

def get_num_of_exited_states(resstring):
    return sum([int(x) for x in resstring])

def get_exited_state_dist(result, num_qubits):
    
    counts = dict(zip(
        [x[:num_qubits] for x in result.get_counts().keys()], 
        [x for x in result.get_counts().values()]
    ))

    cd = {}
    for k,v in counts.items():
        cd[k] = [get_num_of_exited_states(k),v]   

    exi_counts = dict(zip([x for x in range(num_qubits)], [0 for x in range(num_qubits)]))
    for i in range(num_qubits):
        exi_counts[i] = sum([x[1] for x in cd.values() if x[0]==i])
        
    return exi_counts