"""
File name: 
    OTOCProtocol.py
Author: 
    J. Abraham Hernández
Programme function: 
    Outputs the time dependent behaviour of the out-of-time-ordered correlation function as the unitary evolution of a "butterfly" operator is calculated 
    Takes a unitary many-body scrambling operator as input
"""

from qiskit import *
from qiskit import Aer, transpile
import qiskit.quantum_info as qi
from ScramblingUnitary import *
import numpy as np
import matplotlib.pyplot as plt


# Circuit preparation
def ButterflyCircuit(n, k, qubitB, notRef=True):
    """
    Return the resulting "Butterfly" circuit (BC) where the Scrambling unitary was used
    Args:
        n (int): number of qubits
        k (int): number of cycles within BC
        bq (int): position of the "Butterfly" qubit, values from 2-(n-1) possible
    Returns: 
        Circuit object BC
        String containing BDM scores for ScrU
    """    
    # Circuit initialization
    bqr = QuantumRegister(n,"q")
    bc = QuantumCircuit(bqr)
    # Qubit preparation
    bc.h([bqr[i] for i in range(n)])
    bc.s(0)
    bc.cz(1,0)
    # Scrambling unitary initialization
    scramblingO, bdmscore = ScramblingU(n, k, nearest=False, filters=False, printArray=False)
    # ... and its hermitian conjugate
    scramblingOdagg = Udagger(scramblingO)
    # Scrambling unitary, butterfly qubit and scrambling inverse
    bc.append(scramblingO, [bqr[i] for i in range(n)])
    if notRef:
        bc.x(qubitB)
    bc.append(scramblingOdagg, [bqr[i] for i in range(n)])
    # Phase collection / preparation for readout
    # TODO: does it matter in which order the controlled operations are done? i.e. are they symmetric?
    bc.cz(1,0)
    # Not needed as measurement is done natively in the y-basis
    # bc.sdg(0)
    # bc.h(0)
    # BDM score output
    finalBDMscore = "".join("BDM for U: {v1}\n".format(v1=bdmscore))
    
    return bc, finalBDMscore

def OTOCoutcome(circ):
    """ OTOC measurement via expval of Pauli Y on the ancilla qubit
    Args:
        circ (Circuit): valid quantum circuit
    Returns: 
        Float: expectation value of Pauli Y, in [-1,1]
    """    
    # Noisy simulator definition
    simulator = Aer.get_backend("aer_simulator")
    otocCirc = transpile(circ, simulator)
    # Operator to take expectation value of
    qubitM = qi.Pauli("Y")
    qubitT = [0]
    otocCirc.save_expectation_value(qubitM, qubitT)
    result = simulator.run(otocCirc).result()
    exp = result.data()["expectation_value"]
    return exp

def avgOTOC(cyclesK, statRep, qubits):
    """ Average OTOC statistics and measurement results
    Args:
        cyclesK (int): cycle of U, quantifies depth of U
        statRep (int): number of measurement repetition (random unitary generation) 
        qubits (int): total number of qubits in the circuit
    Returns: 
        Plot. Every curve denotes a butterfly position, the x-axis is the depth of U and the y-axis is the expval of Pauli Y
    """    
    # Allowed indices for the butterfly qubit
    qubitBindices = range(1, qubits)
    # Prealocating an array for expectation values
    expvalPauliY = np.empty((qubits, cyclesK), int)
    for n in qubitBindices:
        EVperCycle = []
        for l in range(1, cyclesK+1):
            normVal = []
            for m in range(statRep+1):
                refCirc, refBDM = ButterflyCircuit(qubits, l, n, notRef=False)
                expCirc, expBDM = ButterflyCircuit(qubits, l, n, notRef=True)
                refVal = OTOCoutcome(refCirc)
                expVal = OTOCoutcome(expCirc)
                normVal.append(expVal/refVal)
            normValTotal = np.average(normVal)
            EVperCycle.append(normValTotal)
        expvalPauliY[n-1] = EVperCycle
    # return expvalPauliY

    # Plotting the resulting 2D array
    for i in range(qubits):
        plt.plot(np.arange(1, cyclesK+1, 1, dtype=int), expvalPauliY[i,:], '.-', markersize=12)
    # plt.xlim([1.5, 12.5])
    # plt.ylim([1.5, 4.5])
    plt.legend(['Butterfly correlation decay'])
    plt.xlabel('Number of cycles in $U$', fontdict={'size':20})
    plt.ylabel('$\overline{\text{OTOC}}$', fontdict={'size':20})
    plt.tick_params(axis='x', labelsize=12)
    plt.tick_params(axis='y', labelsize=12)
    plt.show()

if __name__ == "__main__":

    # --- Example implementation ---

    # Fire up the functions
    butterflycircuit, score = ButterflyCircuit(5,1,3, notRef=True)

    # Show the results
    # print(butterflycircuit)
    # print(score)
    # for i in range(5):
    #     butterflycircuit, score = ButterflyCircuit(5,1,3, notRef=True)
    #     print(OTOCoutcome(butterflycircuit))

    # print([j for j in range(1,3+1)])

    # Avg OTOC calculation
    # avgOTOC(2,2,3)
    # print(avgOTOC(3,1,5))
    # print(avgOTOC(3,1,5).shape)
    # arr0 = [i for i in range(1,6)] 
    # arr1 = [i for i in range(2,7)] 
    # arr2 = [i for i in range(1,6)] 
    # arr3 = [i for i in range(2,7)]
    # arr4 = [i for i in range(1,6)]
    # arr5 = np.vstack((arr0,arr1,arr2,arr3,arr4))

    # arr6 = [i for i in range(5)]
    # # arr1 = np.array([i for i in range(1,4)], [j for j in range(1,3)])
    # print(arr5[1,:])
    # print(np.arange(1, 6, 1, dtype=int))
    # # for i in range(5+1):
    # plt.plot(np.arange(1, 6, 1, dtype=int), arr5[1,:], '.-', markersize=12)
    # plt.show()
