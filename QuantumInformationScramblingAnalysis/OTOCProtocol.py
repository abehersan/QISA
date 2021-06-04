"""
File name: 
    OTOCProtocol.py
Author: 
    J. Abraham Hernández
Programme function: 
    Outputs the time dependent behaviour of the out-of-time-ordered correlation function as the unitary evolution of a "butterfly" operator is calculated 
    Takes a unitary many-body scrambling operator as input
"""

from ScramblingUnitary import *
from numpy.core.fromnumeric import shape
from qiskit import *
from qiskit import Aer, transpile
import qiskit.quantum_info as qi
import numpy as np
import matplotlib.pyplot as plt
import timeit
from concurrent.futures import ThreadPoolExecutor

# Circuit preparation
def ButterflyCircuit(n, k, qubitB, notRef=True):
    """
    Return the resulting "Butterfly" circuit (BC) where the Scrambling unitary was used
    Args:
        n (int): number of qubits
        k (int): number of cycles within BC
        qubitB (int): position of the "Butterfly" qubit, values from 2-(n-1) possible
        notRef (bool): if True, it simply applies UU^{dagger}=1
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
    bc.cz(0,1)
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
    bc.cz(0,1)
    # Not needed as measurement is done natively in the y-basis
    # bc.sdg(0)
    # bc.h(0)
    # BDM score output
    # finalBDMscore = "".join("BDM for U: {v1}\n".format(v1=bdmscore))
    
    return bc, bdmscore

def OTOCoutcome(circ):
    """ OTOC measurement via expval of Pauli Y on the ancilla qubit
    Args:
        circ (Circuit): valid quantum circuit
    Returns: 
        Float: expectation value of Pauli Y, in [-1,1]
    """    
    # Noisy simulator definition
    simulator = Aer.get_backend("qasm_simulator")
    otocCirc = transpile(circ, simulator)
    # Operator to take expectation value of
    qubitM = qi.Pauli("Y")
    qubitT = [0]
    otocCirc.save_expectation_value(qubitM, qubitT)
    result = simulator.run(otocCirc).result()
    exp = result.data()["expectation_value"]
    return exp

def avgOTOC(cyclesK, statRep, qubits, graph = True, absVal = True):
    """ Average OTOC statistics and measurement results
    Args:
        cyclesK (int): cycle of U, quantifies depth of U
        statRep (int): number of measurement repetition (random unitary generation) 
        qubits (int): total number of qubits in the circuit
        graph (bool): plotted or not
        absVal (bool): decide if you proceed with OTOC or |OTOC|
    Returns: 
        Plot if graph switch is True. Every curve denotes a butterfly position, the x-axis is the depth of U and the y-axis is the expval of Pauli Y
        Array if graph switch if False, for further processing
    """    
    # Allowed indices for the butterfly qubit
    qubitBindices = range(1, qubits)
    # Prealocating an array for expectation values
    expvalPauliY = np.zeros(shape=(qubits, cyclesK))
    totBDM = np.zeros(shape=(qubits, cyclesK))
    # expvalPauliY = []
    for n in qubitBindices:
        for l in range(1, cyclesK):
            avgBDM = []
            normVal = []
            for m in range(statRep+1):
                # refCirc, refBDM = ButterflyCircuit(qubits, l, n, notRef=False)
                expCirc, expBDM = ButterflyCircuit(qubits, l, n, notRef=True)
                # refVal = OTOCoutcome(refCirc)
                expVal = OTOCoutcome(expCirc)
                # print("buttqub in pos {}, expval {} ".format(n, expVal))
                # normVal.append(expVal/refVal)
                normVal.append(expVal)
                avgBDM.append(expBDM)
            BDMValTotal = np.average(avgBDM)
            normValTotal = np.average(normVal)
            # Absolute value switch
            # if absVal:
            #     expvalPauliY[n][l]=abs(normValTotal) 
            # else:
            expvalPauliY[n][l]=normValTotal 
            totBDM[n][l]=BDMValTotal
        # print(expvalPauliY)
    if absVal:
        expvalPauliY = np.absolute(expvalPauliY)
    print(expvalPauliY.shape)
    if graph:
        # Plotting the resulting 2D arrays
        # General themes
        with plt.style.context("dark_background"):
        # OTOC plot
            ax1 = plt.subplot(211)
            for i in range(qubits):
                plt.plot(np.arange(1, cyclesK+1, 1, dtype=int), expvalPauliY[i,:], '.-', markersize=12)
            plt.setp(ax1.get_xticklabels(), visible=False)
            # Title and text
            plt.title('BUTTERFLY EFFECT IN A QUANTUM CIRCUIT')
            # legends = ["B in qubit {}".format(j) for j in qubitBindices]
            # plt.legend(legends)
            # plt.xlabel('Number of cycles in $U$', fontdict={'size':20})
            plt.ylabel('$|\overline{OTOC}}|$', fontdict={'size':20})
            plt.tick_params(axis='x', labelsize=12)
            plt.tick_params(axis='y', labelsize=12)
            
            # BDM plot
            ax2 = plt.subplot(212, sharex = ax1)
            for i in range(qubits):
                plt.plot(np.arange(1, cyclesK+1, 1, dtype=int), totBDM[i,:], '.-', markersize=12)
            plt.setp(ax2.get_xticklabels(), visible=True, fontsize=8)
            # Title and text
            # plt.title('Butterfly correlation function decay')
            legends = ["B in qubit {}".format(j) for j in qubitBindices]
            plt.legend(legends)
            plt.xlabel('Number of cycles in $U$', fontdict={'size':20})
            plt.ylabel('$\overline{BDM}}$', fontdict={'size':20})
            plt.tick_params(axis='x', labelsize=12)
            plt.tick_params(axis='y', labelsize=12)
            # plt.savefig("fig_1")
            plt.show()
    return expvalPauliY, BDMValTotal
        

if __name__ == "__main__":

    # --- Example implementation ---

    # - Timer starts - 
    t_0 = timeit.default_timer()

    # Fire up a single shot function
    butterflycircuit, score = ButterflyCircuit(5,1,3, notRef=True)

    # Show the results
    print("The OTOC measuring circuit looks like this!\n")
    print(butterflycircuit)
    print("Its characteristic BDM score turned out to be: {}\n".format(score))

    # Avg OTOC calculation for many statistical repetitions
    with ThreadPoolExecutor(4) as ex:
        # ex.map(avgOTOC(51,51,5, True, True))
        print("In a moment, you will effectively see the OTOC(t)'s time-dependence\n")
        print("Circuit depth, statistical repetitions and the number of qubits in our noisy channel contribute to the time it takes for the OTOC to evaluate.\n")
        print("Generally it should take from 1-30 minutes in order to get meaningful (and publishable) results.\n")
        print("The output plot will pop up whenever the computation is done! (:\n")
        ex.map(avgOTOC(20,30,5, True, True)) # produced fig 1
        # ex.map(avgOTOC(20,30,5, True, False)) # produced fig 2

    # - Time checkpoint - 
    t_1 = timeit.default_timer() - t_0
    print("Exec time: {}min".format(t_1/60))
    