from pauli_list import Pauli

from openfermion import fermi_hubbard, hermitian_conjugated
from openfermion.ops import FermionOperator
from openfermion.transforms import jordan_wigner, bravyi_kitaev, bravyi_kitaev_fast
import itertools

def fermiop_Pn(n):
    '''
    Path graph on n vertices, nearest neighbor interactions
    '''
    fermion_operator = FermionOperator()
    for i in range(n-1):
        fermion_operator += FermionOperator('%d^ %d' % (i,i+1), 1.0) + FermionOperator('%d^ %d' % (i+1,i), 1.0)
    return fermion_operator

def fermiop_Cn(n):
    '''
    Cycle graph on n vertices, nearest neighbor interactions
    '''
    return fermiop_Pn(n) + FermionOperator('%d^ %d' % (0,n-1), 1.0) + FermionOperator('%d^ %d' % (n-1, 0), 1.0)

def fermiop_Kn(n):
    '''
    Complete graph on n vertiecs, all-to-all interactions
    '''
    fermion_operator = FermionOperator()
    for i in range(n):
        for j in range(i+1, n):
            fermion_operator += FermionOperator('%d^ %d' %(i,j), 1.0) + FermionOperator('%d^ %d' % (j,i), 1.0)
    return fermion_operator

def single_ops(n):
    fermion_operator = FermionOperator()
    for i in range(n):
        fermion_operator += FermionOperator('%d' % i, 1.0)
    return fermion_operator

def fermiop_line_range(n, max_range, include_self = False):
    fermion_operator = FermionOperator()
    for i in range(n):
        if include_self:
            fermion_operator += FermionOperator('%d^ %d' % (i, i), 1.0)
        for j in range(i+1, min(i+max_range+1, n)):
            fermion_operator += FermionOperator('%d^ %d' %(i,j), 1.0) + FermionOperator('%d^ %d' % (j,i), 1.0)
    return fermion_operator

def fermiop_nn_grid(n_arr, density_interactions=False):
    '''
    n_arr: 2D array with the numbering of the sites
    '''
    fermion_operator = FermionOperator()
    for i in range(len(n_arr)):
        for j in range(len(n_arr[0]) - 1):
            fermion_operator += FermionOperator('%d^ %d' % (n_arr[i][j], n_arr[i][j+1]), 1.0)
            fermion_operator += FermionOperator('%d^ %d' % (n_arr[i][j+1], n_arr[i][j]), 1.0)
            if density_interactions:
                fermion_operator += FermionOperator('%d^ %d^ %d %d' % ((n_arr[i][j], n_arr[i][j+1]) * 2), 1.0)
    for i in range(len(n_arr)-1):
        for j in range(len(n_arr[0])):
            fermion_operator += FermionOperator('%d^ %d' % (n_arr[i][j], n_arr[i+1][j]), 1.0)
            fermion_operator += FermionOperator('%d^ %d' % (n_arr[i+1][j], n_arr[i][j]), 1.0)
            if density_interactions:
                fermion_operator += FermionOperator('%d^ %d^ %d %d' % ((n_arr[i][j], n_arr[i+1][j]) * 2), 1.0)
    return fermion_operator

def fermiop_hubbard_grid(up_arr, down_arr):
    '''
    {up,down}_arr: 2D array with the numbering of the sites, should be same shape
    '''
    length = len(up_arr)
    width = len(up_arr[0])
    assert len(down_arr) == length and len(down_arr[0]) == width
    N = length * width
    fermion_operator = FermionOperator()
    for i in range(length):
        for j in range(width - 1):
            fermion_operator += FermionOperator('%d^ %d' % (up_arr[i][j], up_arr[i][j+1]), 1.0)
            fermion_operator += FermionOperator('%d^ %d' % (up_arr[i][j+1], up_arr[i][j]), 1.0)
            fermion_operator += FermionOperator('%d^ %d' % (down_arr[i][j], down_arr[i][j+1]), 1.0)
            fermion_operator += FermionOperator('%d^ %d' % (down_arr[i][j+1], down_arr[i][j]), 1.0)
    for i in range(length-1):
        for j in range(width):
            fermion_operator += FermionOperator('%d^ %d' % (up_arr[i][j], up_arr[i+1][j]), 1.0)
            fermion_operator += FermionOperator('%d^ %d' % (up_arr[i+1][j], up_arr[i][j]), 1.0)
            fermion_operator += FermionOperator('%d^ %d' % (down_arr[i][j], down_arr[i+1][j]), 1.0)
            fermion_operator += FermionOperator('%d^ %d' % (down_arr[i+1][j], down_arr[i][j]), 1.0)
    for i in range(length):
        for j in range(width):
            fermion_operator += FermionOperator('%d^ %d^ %d %d' % ((up_arr[i][j], down_arr[i][j]) * 2), 1.0)
    return fermion_operator

################# OLDER

def fermiop_hubbard(n_sites):
    '''
    Hubbard model on n sites, restricted lattice structure
    TODO describe
    TODO paramters
    '''
    U = 2.0
    J = -1.0
    hubbard = fermi_hubbard(1, n_sites, tunneling=-J, coulomb=U, periodic=False)
    return hubbard

def fermiop_hubbard_all(n_sites):
    '''
    TODO
    '''
    n_qubit = n_sites*2
    ops = [FermionOperator(((i, 1), (j, 0)), coefficient =1) + hermitian_conjugated(FermionOperator(((i, 1), (j, 0)), coefficient =1)) for i in range(n_qubit - 2) for j in range(i+2,n_qubit, 2)]
    ham=sum(ops)
    return ham

def fermiop_hubbard_all_spinless(n_sites):
    n_qubit = n_sites
    ops = [FermionOperator(((i, 1), (j, 0)), coefficient =1) + hermitian_conjugated(FermionOperator(((i, 1), (j, 0)), coefficient =1)) for i in range(n_qubit - 1) for j in range(i+1,n_qubit, 1)]
    ham=sum(ops)
    return ham


def fermiop_hubbard_nn_all(n_sites): #Plus density-density interaction
    n_qubit = n_sites*2
    ops = [FermionOperator(((i, 1), (j, 0)), coefficient =1) + hermitian_conjugated(FermionOperator(((i, 1), (j, 0)), coefficient =1)) + FermionOperator(((i, 1), (i, 0),(j, 1), (j, 0)), coefficient =0.1212) + hermitian_conjugated(FermionOperator(((i, 1), (i, 0),(j, 1), (j, 0)), coefficient =0.1212)) for i in range(n_qubit - 2) for j in range(i+2,n_qubit, 2)]
    ham=sum(ops)
    return ham

def fermiop_hubbard_aaaa_all_old(n_sites): #Plus all two-body interaction
    n_qubit = n_sites*2
    ops = [FermionOperator(((i, 1), (j, 0)), coefficient=1) + hermitian_conjugated(FermionOperator(((i, 1), (j, 0)), coefficient=1)) for i in range(n_qubit - 2) for j in range(i+2, n_qubit, 2)]
    ops += [FermionOperator(((i, 1), (j, 0), (k, 1), (l, 0)), coefficient=0.1212) + hermitian_conjugated(FermionOperator(((i, 1), (j, 0), (k, 1), (l, 0)), coefficient=0.1212)) for i in range(n_qubit - 4) for j in range(i+2, n_qubit, 2) for k in range(j+2, n_qubit, 2)for l in range(k+2, n_qubit, 2)]
    ham=sum(ops)
    #print(ham)
    return ham

def fermiop_hubbard_aaaa_all(n_sites): #Plus all two-body interaction
    n_qubit = n_sites
    ops = [FermionOperator(((i, 1), (j, 0)), coefficient=1) + hermitian_conjugated(FermionOperator(((i, 1), (j, 0)), coefficient=1)) for i in range(n_qubit - 1) for j in range(i+1, n_qubit,1)]
    #ops += [FermionOperator(((i, 1), (j, 0), (k, 1), (l, 0)), coefficient=0.1212) + hermitian_conjugated(FermionOperator(((i, 1), (j, 0), (k, 1), (l, 0)), coefficient=0.13)) for i in range(n_qubit - 1) for j in range(n_qubit-1) for k in range(i+1,n_qubit,1)for l in range(j+1, n_qubit, 1)]
    ops += [FermionOperator(((i, 1), (j, 0), (k, 1), (l, 0)), coefficient=0.1212)  for i in range(n_qubit - 1) for j in range(n_qubit-1) for k in range(i+1,n_qubit,1)for l in range(j+1, n_qubit, 1)]
    ham=sum(ops)
    print(ham)
    return ham


def fermiop_nn_triangular(n_sites): # a_i^d a_j + n_i m_j for a triangular lattice. n_site should be a multiple of xaxis
    n_qubit = n_sites*2
    xaxis=4
    yaxis=int(n_sites/xaxis)
    all_vertices = list(itertools.product(range(xaxis),range(yaxis)))
    print(all_vertices)
    ops=FermionOperator(term=None, coefficient=1.0)
    for vertex in all_vertices:
        print(vertex, vertex[0], vertex[1] + 1)
        # Neighboring vertices
        vertex_n = [(vertex[0]+1,vertex[1]),(vertex[0]-1,vertex[1]+1), (vertex[0],vertex[1]+1)]
        for i in range(3):
            vertex_n1=vertex_n[i]
            if vertex_n1 in all_vertices:
                print("vertex_n1",vertex_n1)
                ops += FermionOperator(((vertex[0]+vertex[1]*xaxis, 1), (vertex_n1[0]+vertex_n1[1]*xaxis, 0)), coefficient =1) + hermitian_conjugated(FermionOperator(((vertex[0]+vertex[1]*xaxis, 1), (vertex_n1[0]+vertex_n1[1]*xaxis, 0)), coefficient =1))
    return ops


def qubitop_to_pauli(qubit_op, num_qubits):
    '''
    Convert openfermion.ops.operator.qubit_operator.QubitOperator to Pauli
    TODO: infer num_qubits
    '''
    pauli = Pauli(num_qubits)
    terms = qubit_op.terms
    for term in terms:
        pauli_string = [0]*num_qubits
        for (i, c) in term:
            pauli_string[i] = 1 if c=='X' else 2 if c=='Y' else 3 if c=='Z' else "ERROR"
        pauli.add(tuple(pauli_string), terms[term])
    return pauli

def fermiop_to_pauli(fermi_op, num_qubits, transform=jordan_wigner):
    '''
    Convert openfermion.ops.FermionOperator to Pauli using qubit transform
    TODO: infer num_qubits
    '''
    return qubitop_to_pauli(transform(fermi_op), num_qubits)
