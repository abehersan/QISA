import numpy as np
import qiskit as q

from qiskit.pulse import *

def update_basis_gates_and_circ_inst_map(decomposed_circuit, backend, circ_inst_map):
    """While Parametrized Schedules are not supported, simply update basis gates and circ_inst_map."""
    basis_gates = backend.configuration().basis_gates
    for instruction, qargs, cargs in decomposed_circuit.data:
        if instruction.name.startswith('direct_rx'):
            if instruction.name not in basis_gates:
                basis_gates.append(instruction.name)
            
            if not circ_inst_map.has(instruction.name, qubits=[qargs[0].index]):
                theta = float(instruction.name[len('direct_rx_'):])
                schedule = get_direct_rx_schedule(theta, qargs[0].index, circ_inst_map, backend)
                circ_inst_map.add(instruction.name, qubits=[qargs[0].index], schedule=schedule)

        elif instruction.name.startswith('cr'):
            if instruction.name not in basis_gates:
                basis_gates.append(instruction.name)

            if not circ_inst_map.has(instruction.name, qubits=[qargs[0].index, qargs[1].index]):
                theta = float(instruction.name[len('cr_'):])
                schedule = get_cr_schedule(theta, qargs[0].index, qargs[1].index, circ_inst_map, backend)
                circ_inst_map.add(instruction.name,
                            qubits=[qargs[0].index, qargs[1].index], schedule=schedule)

        elif instruction.name == 'open_cx':
            if instruction.name not in basis_gates:
                basis_gates.append(instruction.name)


def rescale_samples(samples, scale_factor, method='rescale_height'):
    assert scale_factor <= 1, 'only tested for scaling down pulses'

    if method == 'rescale_height':
        return _rescale_height(samples, scale_factor)
    elif method == 'rescale_width':
        return _rescale_width(samples, scale_factor)
    elif method == 'rescale_height_and_width':
        return _rescale_height_and_width(samples, scale_factor)


def _rescale_height(samples, scale_factor):
    return samples * scale_factor


def _rescale_width(samples, scale_factor):
    assert False, 'still debugging implementation'


def _rescale_height_and_width(samples, scale_factor):
    assert False, 'still debugging implementation'

    print('original real area under curve is %s' % sum(map(np.real, samples)))

    rescaled_length = int(0.5 + len(samples) * np.sqrt(scale_factor))
    rescaled_samples = [0] * rescaled_length
    width_scale_factor = rescaled_length / len(samples)
    height_scale_factor = scale_factor / width_scale_factor
    samples = samples * scale_factor
    
    for i in range(len(samples)):
        # split samples[i...i+1] into rescaled_samples[i/scale_factor...(i+1)/scale_factor]
        if int(i * width_scale_factor) == int((i + 1) * width_scale_factor):
            rescaled_samples[int(i * width_scale_factor)] += samples[i]
        else:
            fraction = int(1 + i * width_scale_factor) - int(i * width_scale_factor)
            rescaled_samples[int(i * width_scale_factor)] += samples[i] * fraction
            rescaled_samples[int(i * width_scale_factor)] += samples[i] * (1 - fraction)
    print('final real area under curve is %s' % sum(map(np.real, rescaled_samples)))
    return rescaled_samples


def get_direct_rx_schedule(theta, qubit_index, backend):
    circ_inst_map = backend.defaults().instruction_schedule_map 
    
    x_instructions = circ_inst_map.get('x', qubits=[qubit_index]).instructions
    assert len(x_instructions) == 1
    x_samples = x_instructions[0][1].command.samples
    area_under_curve = sum(map(np.real, x_samples))
    
    if theta > np.pi:
        theta -= 2 * np.pi
    
    direct_rx_samples = rescale_samples(x_samples, (theta / np.pi))
    direct_rx_samplepulse = q.pulse.SamplePulse(direct_rx_samples)
    direct_rx_command = direct_rx_samplepulse(backend.configuration().drive(qubit_index))
    return q.pulse.Schedule([0, direct_rx_command])

    return direct_rx_command
