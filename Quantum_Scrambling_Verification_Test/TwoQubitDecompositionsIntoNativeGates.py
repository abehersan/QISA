#  Source:
#  "Optimized Quantum Compilation for Near-Term Algorithms with OpenPulse"; arXiv:2004.11205 [quant-ph]; 
#  Source code:
#  https://github.com/singular-value/optimizations_via_openpulse/blob/master/TwoQubitDecompositionsIntoNativeGates.ipynb

import numpy as np
from scipy.optimize import minimize, NonlinearConstraint
from scipy.linalg import fractional_matrix_power, expm
from qiskit.quantum_info.synthesis.two_qubit_decompose import TwoQubitBasisDecomposer
import qiskit as q

CNOT_UNITARY = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
ISWAP_UNITARY = np.array([[1, 0, 0, 0], [0, 0, -1j, 0], [0, -1j, 0, 0], [0, 0, 0, 1]])
CR90_UNITARY = np.array([[1, -1j, 0, 0], [-1j, 1, 0, 0], [0, 0, 1, 1j], [0, 0, 1j, 1]]) / np.sqrt(2)

class NativeTwoQubitGate(object):
    num_params = -1

    def get_unitary(self, params=[]):
        assert len(params) == self.num_params
        assert self.num_params != -1, 'Subclass should set num_params'
        assert self.num_params in [0, 1], 'Currently doesn\'t handle multi-parameter native gates'

        return self._get_unitary(params)
    
    def _get_unitary(self, params=[]):
        raise NotImplemented('Subclass should implement')
        
    def __str__(self):
        return self.__class__.__name__

    
class NativeCNOT(NativeTwoQubitGate):
    num_params = 0
    
    def _get_unitary(self, params=[]):
        return CNOT_UNITARY


class NativeiSWAP(NativeTwoQubitGate):
    num_params = 0
    
    def _get_unitary(self, params=[]):
        return ISWAP_UNITARY
    

class NativeRootiSWAP(NativeTwoQubitGate):
    num_params = 0
    
    def _get_unitary(self, params=[]):
        return fractional_matrix_power(ISWAP_UNITARY, 0.5)



class NativeCR(NativeTwoQubitGate):
    num_params = 1
    
    def _get_unitary(self, params=[]):
        theta = params[0]
        return np.array([
            [np.cos(theta/2), -1j * np.sin(theta/2), 0, 0],
            [-1j * np.sin(theta/2), np.cos(theta/2), 0, 0],
            [0, 0, np.cos(theta/2), 1j * np.sin(theta/2)],
            [0, 0, 1j * np.sin(theta/2), np.cos(theta/2)]])


class NativeCZ(NativeTwoQubitGate):
    num_params = 0
    
    def _get_unitary(self, params=[]):
        return np.diag([1,1,1,-1])
    
    
class NativeCR90(NativeTwoQubitGate):
    num_params = 0
    
    def _get_unitary(self, params=[]):
        return CR90_UNITARY


class NativeParametrizediSWAP(NativeTwoQubitGate):
    num_params = 1
    
    def _get_unitary(self, params=[]):
        theta = params[0]
        return np.array([
            [1, 0, 0, 0],
            [0, np.cos(theta), -1j*np.sin(theta), 0],
            [0, -1j*np.sin(theta), np.cos(theta), 0],
            [0, 0, 0, 1]])
    
    
class NativeSWAPAlpha(NativeTwoQubitGate):
    num_params = 1
    
    def _get_unitary(self, params=[]):
        alpha = params[0]
        return np.array([
            [1, 0, 0, 0],
            [0, (1 + np.exp(1j*np.pi*alpha))/2, (1 - np.exp(1j*np.pi*alpha))/2, 0],
            [0, (1 - np.exp(1j*np.pi*alpha))/2, (1 + np.exp(1j*np.pi*alpha))/2, 0],
            [0, 0, 0, 1]
        ])

    
class NativeBSWAP(NativeTwoQubitGate):
    num_params = 0
    
    def _get_unitary(self, params=[]):
        return np.array([
            [np.cos(np.pi/8), 0, 0, 1j * np.sin(np.pi/8)],
            [0, np.cos(3*np.pi/8), 1j * np.sin(3*np.pi/8), 0],
            [0, 1j * np.sin(3*np.pi/8), np.cos(3*np.pi/8), 0],
            [1j * np.sin(np.pi/8), 0, 0, np.cos(np.pi/8)]
        ])


class NativeMAP(NativeTwoQubitGate):
    num_params = 0

    def _get_unitary(self, params=[]):
        return expm(-1j * np.pi/4 * np.diag([1, -1, -1, 1]))

    
class NativeRIP(NativeTwoQubitGate):
    pass


class NativeMS(NativeTwoQubitGate):
    num_params = 1
    
    def _get_unitary(self, params=[]):
        theta = params[0]
        return np.array([
            [np.cos(np.pi*theta/2), 0, 0, -1j * np.sin(np.pi*theta/2)],
            [0, np.cos(np.pi*theta/2), -1j*np.sin(np.pi*theta/2), 0],
            [0, -1j*np.sin(np.pi*theta/2), np.cos(np.pi*theta/2), 0],
            [-1j*np.sin(np.pi*theta/2), 0, 0, np.cos(np.pi*theta/2)]
        ])
    
class TargetTwoQubitOperation(object):
    def get_unitaries(self, params=None):
        raise NotImplemented('Subclass should implement')
        
    def __str__(self):
        return self.__class__.__name__


class TargetCNOT(TargetTwoQubitOperation):
    def get_unitaries(self):
        return [CNOT_UNITARY]
    

class TargetSWAP(TargetTwoQubitOperation):
    def get_unitaries(self):
        return [np.array([[1, 0, 0, 0],
                          [0, 0, 1, 0],
                          [0, 1, 0 ,0],
                          [0, 0, 0, 1]])]

    
class TargetZZInteraction(TargetTwoQubitOperation):
    def get_unitaries(self, params=None):
        if params is None:
            params = np.random.random(10) * 2 * np.pi
        return [np.diag([1, np.exp(1j * param), np.exp(1j * param), 1]) for param in params]
    
    
class TargetZZSWAP(TargetTwoQubitOperation):
    def get_unitaries(self, params=None):
        if params is None:
            params = np.random.random(10) * 2 * np.pi
        unitaries = []
        for param in params:
            unitaries.append(np.array([[0, 0, 0, 1],
                                       [0, np.exp(1j*param), 0, 0],
                                       [0, 0, np.exp(1j*param), 0],
                                       [1, 0, 0, 0]]))
        return unitaries
    

class TargetFermionicSimulation(TargetTwoQubitOperation):
    def get_unitaries(self, params=None):
        params = [np.random.random(2) * 2 * np.pi for _ in range(20)]
        unitaries = []
        for T, V in params:
            unitaries.append(np.array([[1, 0, 0, 0],
                                       [0, -1j*np.sin(T), np.cos(T), 0],
                                       [0, np.cos(T), -1j*np.sin(T), 0],
                                       [0, 0, 0, -np.exp(-1j*V)]]))
        return unitaries
    
    
class TargetFermionicFourierTransform(TargetTwoQubitOperation):
    def get_unitaries(self):
        return [np.array([
            [1, 0, 0, 0],
            [0, 1/np.sqrt(2), 1/np.sqrt(2), 0],
            [0, 1/np.sqrt(2), -1/np.sqrt(2), 0],
            [0, 0, 0, -1]
        ])]
    

class TargetBogoliubovTransform(TargetTwoQubitOperation):
    def get_unitaries(self, params=None):
        if params is None:
            params = np.random.random(10) * 2 * np.pi

        def _get_unitary(expo):  # adapted from Cirq
            # --X--S--|iSWAP^expo|--S^1.5--X--
            # --------|iSWAP^expo|------------
            U = np.kron(np.array([[0, 1], [1, 0]]), np.eye(2))
            U = np.kron(fractional_matrix_power(np.diag([1, -1]), 1.5), np.eye(2)) @ U
            U = fractional_matrix_power(ISWAP_UNITARY, expo) @ U
            U = np.kron(np.array(np.diag([1, -1])), np.eye(2)) @ U
            U = np.kron(np.array([[0, 1], [1, 0]]), np.eye(2)) @ U
            return U

        return [_get_unitary(param) for param in params]