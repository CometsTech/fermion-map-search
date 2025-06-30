import numpy as np
import scipy.optimize
import os

PAULI_MATRICES = (np.eye(2), np.array([[0, 1], [1,0]]), np.array([[0, (0-1j)], [(0+1j) ,0]]), np.array([[1, 0], [0, -1]]))

class Pauli:
    def __init__(self, num_qubits, terms=None, prec=1e-6):
        self.num_qubits = num_qubits
        # term = tuple of numbers representing tensor product, 0=I, 1=X, 2=Y, 3=Z
        self.terms = terms or {} # term : coef
        self.length = len(self.terms)
        self.prec = prec

    def add(self, term, coef, unique=False):
        if not unique and term in self.terms:
            self.terms[term] += coef
            if abs(self.terms[term]) < self.prec:
                self.terms.pop(term)
                self.length -= 1
        else:
            self.terms[term] = coef
            self.length += 1

    def conj(self):
        # Returns conjugate copy of self
        return Pauli(self.num_qubits, {term: np.conj(coef) for term,coef in self.terms.items()})

    def to_matrix(self):
        # Returns np.array matrix form
        dim = 2**self.num_qubits
        mat = np.zeros((dim, dim), 'complex128')
        for term in self.terms:
            mat_term = np.eye(1)
            for i in term:
                mat_term = np.kron(mat_term, PAULI_MATRICES[i])
            mat += self.terms[term] * mat_term
        return mat

    def __len__(self):
        return self.length
    
    def __str__(self):
        return ''.join([str(np.round(coef, 6)) + ' ' + str(term) + '\n'  for (term, coef) in self.terms.items()])

    def __eq__(self, other):
        # assumes same prec
        if type(other) != Pauli or len(other) != self.length:
            return False
        for term in self.terms:
            if term not in other.terms or abs(self.terms[term] - other.terms[term]) > self.prec:
                print(term)
                return False
        return True

def mult_pauli_terms(pauli1, pauli2):
    # returns (coef, term)
    num_qubits = len(pauli1)
    coef = 1.0
    term = [0]*num_qubits
    for i in range(num_qubits):
        if pauli1[i] == 0:
            term[i] = pauli2[i]
        elif pauli2[i] == 0:
            term[i] = pauli1[i]
        elif pauli1[i] == pauli2[i]:
            term[i] = 0
        elif pauli1[i] == 1 and pauli2[i] == 2:
            coef *= 1j
            term[i] = 3
        elif pauli1[i] == 2 and pauli2[i] == 3:
            coef *= 1j
            term[i] = 1
        elif pauli1[i] == 3 and pauli2[i] == 1:
            coef *= 1j
            term[i] = 2
        elif pauli1[i] == 3 and pauli2[i] == 2:
            coef *= -1j
            term[i] = 1
        elif pauli1[i] == 2 and pauli2[i] == 1:
            coef *= -1j
            term[i] = 3
        elif pauli1[i] == 1 and pauli2[i] == 3:
            coef *= -1j
            term[i] = 2
        else:
            raise ValueError("HELP pauli_mult")
    return (coef, tuple(term))

def mult_paulis(pauli1, pauli2):
    assert(pauli1.num_qubits == pauli2.num_qubits)
    num_qubits = pauli1.num_qubits
    res = Pauli(num_qubits)
    for (term1, coef1) in pauli1.terms.items():
        for (term2, coef2) in pauli2.terms.items():
            coef, term = mult_pauli_terms(term1, term2)
            res.add(term, coef*coef1*coef2)
    return res

def CNOT(site_C, site_X, num_qubits):
    # Applies a CNOT with control site_C and target site_X
    pauli = Pauli(num_qubits)
    term = [0]*num_qubits
    pauli.add(tuple(term), 0.5)
    term[site_C] = 3
    pauli.add(tuple(term), 0.5)
    term[site_X] = 1
    pauli.add(tuple(term), -0.5)
    term[site_C] = 0
    pauli.add(tuple(term), 0.5)
    return pauli

def H_gate(site, num_qubits):
    pauli = Pauli(num_qubits)
    term = [0]*num_qubits
    term[site] = 1
    pauli.add(tuple(term), 1/np.sqrt(2))
    term[site] = 3
    pauli.add(tuple(term), 1/np.sqrt(2))
    return pauli

def S_gate(site, num_qubits):
    pauli = Pauli(num_qubits)
    term = [0]*num_qubits
    term[site] = 0
    pauli.add(tuple(term), 0.5+0.5j)
    term[site] = 3
    pauli.add(tuple(term), 0.5-0.5j)
    return pauli

def T_gate(site, num_qubits):
    pauli = Pauli(num_qubits)
    term = [0]*num_qubits
    term[site] = 0
    pauli.add(tuple(term), 0.5*(1 + (1 + 1j)/np.sqrt(2)))
    term[site] = 3
    pauli.add(tuple(term), 0.5*(1 - (1+1j)/np.sqrt(2)))
    return pauli

def cost(pauli, eps=1e-6, func=lambda x:1):
    total = 0
    for term, coef in pauli.terms.items():
        if abs(coef) > eps:
            total += func(coef) * np.count_nonzero(term)
    return total

def UAUd(U, A):
    return mult_paulis(mult_paulis(U, A), U.conj())

def UAUd_list(gate_list, H, debug=False):
    num_qubits = H.num_qubits
    A = H
    for gate in gate_list:
        A = UAUd(gate, A)
    if debug:
        print(A, cost(A))
    return A

def UAUd_CNOT_list(pair_list, H, debug=False):
    '''
    Performs CNOT for every pair of qubits in pair_list
    pair_list: list of (control, target) pairs on which to perform CNOT
    '''
    num_qubits = H.num_qubits
    A = H
    for (control, target) in pair_list:
        A = UAUd(CNOT(control, target, num_qubits), A)
    if debug:
        print(A)
        print("cost:", cost(A))
    return A

