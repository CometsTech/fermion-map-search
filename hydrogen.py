from qiskit_nature import settings
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.formats.molecule_info import MoleculeInfo
from qiskit_nature.second_q.mappers import JordanWignerMapper, BravyiKitaevMapper, ParityMapper, BravyiKitaevSuperFastMapper
from qiskit_nature.units import DistanceUnit

settings.use_pauli_sum_op = False

import numpy as np

from pauli_list import *

def hydrogen_line(n, mapper=None, taper=False, debug=False):
    '''
    Returns the qubit operator for second quant for Moluecule with n Hydrogen in chain
    n: number of Hydrogen atoms
    mapper: JordanWignerMapper(), BravyiKitaevMapper(), or ParityMapper()
    if mapper == None: return second_q_op
    '''
    molecule = MoleculeInfo(
        symbols = ['H'] * n,
        coords = [(0, 0, 0.735*i) for i in range(n)],
        multiplicity = 1 + n%2,
        charge = 0,
        units = DistanceUnit.ANGSTROM
    )

    driver = PySCFDriver.from_molecule(molecule=molecule, basis='sto3g')
    problem = driver.run()

    second_q_op = problem.hamiltonian.second_q_op()
    if not mapper:
        return second_q_op
    if taper:
        mapper = problem.get_tapered_mapper(mapper)
    qubit_op = mapper.map(second_q_op)
    if debug:
        return (problem, second_q_op, qubit_op)
    return qubit_op

def op_to_pauli(qubit_op, prec=1e-6):
    # Returns pauli_list.Pauli corresponding to SparsePauliOp operator
    pauli = Pauli(qubit_op.num_qubits, prec=prec)
    for i in range(len(qubit_op.paulis)):
        pauli_string = [0 if c=='I' else 1 if c=='X' else 2 if c=='Y' else 3 if c=='Z' else "ERROR" for c in str(qubit_op.paulis[i])]
        pauli.add(tuple(pauli_string), qubit_op.coeffs[i])
    return pauli
