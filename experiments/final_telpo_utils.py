########################################
# circuit structuring funcs
########################################

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
    
    upper_pairs = [(max(upper_qs), min(upper_qs))] + sorted([(upper_qs[i],upper_qs[i+1]) 
                   for i in range(len(upper_qs)-1)] + [(upper_qs[i],upper_qs[i+2]) 
                                                       for i in range(len(upper_qs)-2)])
    
    lower_pairs = sorted([(lower_qs[i],lower_qs[i+1]) 
                   for i in range(len(lower_qs)-1)] + [(lower_qs[i],lower_qs[i+2]) 
                                                       for i in range(len(lower_qs)-2)])
    lower_pairs = [(min(lower_qs),max(lower_qs))] + [x for x in reversed(lower_pairs)] 

    return upper_pairs, lower_pairs

########################################
# circuit operation funcs
########################################

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