################################################################################
# structuring funcs
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
    
    upper_pairs = [(max(upper_qs), min(upper_qs))] + sorted([(upper_qs[i],upper_qs[i+1]) 
                   for i in range(len(upper_qs)-1)] + [(upper_qs[i],upper_qs[i+2]) 
                                                       for i in range(len(upper_qs)-2)])
    
    lower_pairs = sorted([(lower_qs[i],lower_qs[i+1]) 
                   for i in range(len(lower_qs)-1)] + [(lower_qs[i],lower_qs[i+2]) 
                                                       for i in range(len(lower_qs)-2)])
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
# circuit operation funcs
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
# funcs for general telpo protocol to test scrambling properties
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
        
        
def apply_unitary_test(circ, unitary):
    print("unitary method coming soon")
    print('\n\tsorry, not implemented yet')