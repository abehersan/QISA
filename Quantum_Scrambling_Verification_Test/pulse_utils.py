import numpy as np
from numpy import pi
import qiskit as q
from qiskit.pulse import *

################################################################################
# FUNCS TO GET CR PULSES, PARAMETERS AND UCHANS FROM BACKEND
################################################################################

# get cr pulses from cx basis gates
def get_CRs(qc: int, qt: int, backend):
    inst_sched_map = backend.defaults().instruction_schedule_map 
    cx = inst_sched_map.get('cx', qubits=[qc, qt])
    idx = 0
    # look for first Play instruction on a ControlChannel
    while (type(cx.instructions[idx][1].channels[0]) is not ControlChannel) or \
        (type(cx.instructions[idx][1]) is not Play):
        idx += 1
    return (cx.instructions[idx][1].channels[0], cx.instructions[idx][1].pulse)

# get dict with all cr pulses and params (val) for all coupled pairs (keys)
def get_CR_dict(backend):
    backend_defaults = backend.defaults()
    inst_sched_map = backend_defaults.instruction_schedule_map 
    inst_sched_map.instructions
    backend_config = backend.configuration()
    coupling_map = backend.configuration().coupling_map
    # save control channel nums and CR pulses for all pairs in coupling map to dict
    CR_dict = dict(zip([tuple(pair) for pair in coupling_map], 
             [[get_CRs(pair[0], pair[1],backend)] for pair in coupling_map]
    ))
    return CR_dict

################################################################################
# FUNCS TO GET CR(THETA) AND RX(THETA)
################################################################################

# adapted version of source: 
# "Optimized Quantum Compilation for Near-Term Algorithms with OpenPulse"; arXiv:2004.11205 [quant-ph]; 
# Source code: https://github.com/singular-value/optimizations_via_openpulse

# func to get cr(theta) pulse for qubit pair and backend
def get_cr_schedule(theta, control, target, backend):

    """Returns schedule for a cross-resonance pulse between control and target.
    Does a RX(-theta) on target if control is |0> and a RX(theta) on target if
    control is |1>.
    Crashes if the backend does not support CR between control and target
    (either because no connectivity, or because the CR is between target and control)
    """
    inst_sched_map = backend.defaults().instruction_schedule_map 
    cx_instructions = inst_sched_map.get('cx', qubits=[control, target]).instructions
    xc_instructions = inst_sched_map.get('cx', qubits=[target, control]).instructions
    
    assert len(cx_instructions) < len(xc_instructions), 'CR pulse is on flipped indices'
    
    cr_control_inst, _ = get_CRs(control, target, backend)
    #cr_control_inst = [y for (x,y) in cx_instructions if y.name != None and 'CR90p' in y.name and y.channels[0].name.startswith('u')]
    cr_drive_inst = [y for (x,y) in cx_instructions if y.name != None and 'CR90p' in y.name and y.channels[0].name.startswith('d')]

    #assert len(cr_drive_inst) == 1 and len(cr_control_inst) == 1
    cr_control_inst = cr_control_inst[0]  # driving of control qubit at target's frequency
    cr_drive_inst = cr_drive_inst[0]  # active cancellation tone

    flip = False
    if theta < 0:
        flip = True
        theta = -1 * theta
    
    if theta > 2 * np.pi:
        theta -= 2 * np.pi

    full_area_under_curve = sum(map(np.real, cr_control_inst.pulse.samples))
    target_area_under_curve = full_area_under_curve * (theta / (np.pi / 2))

    # CR pulse samples have gaussian rise, flattop, and then gaussian fall.
    # we want to find the start and end indices of the flattop
    flat_start = 0
    while cr_drive_inst.pulse.samples[flat_start] != cr_drive_inst.pulse.samples[flat_start + 1]:
        flat_start += 1
    assert cr_control_inst.pulse.samples[flat_start] == cr_control_inst.pulse.samples[flat_start + 1]

    flat_end = flat_start + 1
    while cr_drive_inst.pulse.samples[flat_end] == cr_drive_inst.pulse.samples[flat_end + 1]:
        flat_end += 1
    assert cr_control_inst.pulse.samples[flat_end] == cr_control_inst.pulse.samples[flat_end - 1]

    area_under_curve = sum(map(np.real, cr_control_inst.pulse.samples[:flat_start]))
    area_under_curve += sum(map(np.real, cr_control_inst.pulse.samples[flat_end+1:]))
    flat_duration = (target_area_under_curve - area_under_curve) / np.real(cr_control_inst.pulse.samples[flat_start])
    flat_duration = max(0, int(flat_duration + 0.5))
    duration = len(cr_drive_inst.pulse.samples[:flat_start]) + flat_duration + len(cr_drive_inst.pulse.samples[flat_end+1:])
    if duration % 16 <= 8 and flat_duration > 8:
        flat_duration -= duration % 16
    else:
        flat_duration += 16 - (duration % 16)

    cr_drive_samples = np.concatenate([
        cr_drive_inst.pulse.samples[:flat_start],
        [cr_drive_inst.pulse.samples[flat_start]] * flat_duration,
        cr_drive_inst.pulse.samples[flat_end+1:]
    ])

    cr_control_samples = np.concatenate([
        cr_control_inst.pulse.samples[:flat_start],
        [cr_control_inst.pulse.samples[flat_start]] * flat_duration,
        cr_control_inst.pulse.samples[flat_end+1:]
    ])

    assert len(cr_drive_samples) % 16 == 0
    assert len(cr_control_samples) % 16 == 0

    current_area_under_curve = sum(map(np.real, cr_control_samples))
    scaling_factor = target_area_under_curve / current_area_under_curve

    cr_drive_samples *= scaling_factor
    cr_control_samples *= scaling_factor
    
    cr_p_schedule = Play(Waveform(cr_drive_samples),
                         cr_drive_inst.channels[0]) | Play(Waveform(cr_control_samples), cr_control_inst.channels[0])
    cr_m_schedule = Play(Waveform(-1*cr_drive_samples),
                         cr_drive_inst.channels[0]) | Play(Waveform(-1*cr_control_samples),cr_control_inst.channels[0])

    if flip:
        schedule = cr_m_schedule
        schedule |= inst_sched_map.get('x', qubits=[control]) << schedule.duration
        schedule |= cr_p_schedule << schedule.duration
    else:
        schedule = cr_p_schedule
        schedule |= inst_sched_map.get('x', qubits=[control]) << schedule.duration
        schedule |= cr_m_schedule << schedule.duration

    return schedule
