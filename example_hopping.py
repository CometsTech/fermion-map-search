from pauli_list import *
#from hydrogen import *
from models import *
from search import *

from qiskit.quantum_info import SparsePauliOp

NUM_QUBITS = 20
fermion_operator = fermiop_line_range(NUM_QUBITS, 7, False)

bk_operator = bravyi_kitaev(fermion_operator)
jw_operator = jordan_wigner(fermion_operator)
op = qubitop_to_pauli(bk_operator, NUM_QUBITS)
print("Start", cost(op))
#print(op)
#mat = op.to_matrix() # seems hard for n >= 6
#print(mat)
#print(np.linalg.eig(mat)[0])

print(NUM_QUBITS, cost(op), len(op))

ternary_tree_ops = qubitop_to_pauli(bravyi_kitaev(single_ops(NUM_QUBITS)), NUM_QUBITS)
#print(ternary_tree_ops)

#print(search_dfs(op, debug=True))
all_pairs = list(itertools.permutations(range(NUM_QUBITS), 2))
gate_set_1 = {"CNOT(%d, %d)" % (i,j): CNOT(i, j, NUM_QUBITS) for (i, j) in all_pairs}
gate_set_2 = gate_set_1 | {"H(%d) CNOT(%d, %d)" % (i, i, j): mult_paulis(CNOT(i, j, NUM_QUBITS), H_gate(i, NUM_QUBITS)) for (i, j) in all_pairs}
gate_set_2a = gate_set_2 | {"H(%d) CNOT(%d, %d)" % (j, i, j): mult_paulis(CNOT(i, j, NUM_QUBITS), H_gate(j, NUM_QUBITS)) for (i, j) in all_pairs}
gate_set_3 = gate_set_2 | {"S(%d) CNOT(%d, %d)" % (i, i, j): mult_paulis(CNOT(i, j, NUM_QUBITS), S_gate(i, NUM_QUBITS)) for (i, j) in all_pairs}
gate_set_4 = gate_set_1 | {"S(%d) CNOT(%d, %d)" % (i, i, j): mult_paulis(CNOT(i, j, NUM_QUBITS), S_gate(i, NUM_QUBITS)) for (i, j) in all_pairs}
gate_set = gate_set_2
#res = search_astar(op, gate_set, debug=1, timeout=100)
res = search_SA(op, gate_set, lambda n : np.log(1+n/2)/cost(op)*40 if n > 10 else 0, debug=1, max_depth=1_000_000, store_gates=True)
#print(res[0], len(res[1]))
#print(res[2])
#print(UAUd_list([gate_set[i] for i in res[1]]

new_ternary_tree_ops = UAUd_list([gate_set[i] for i in res[1]], ternary_tree_ops)
#print(new_ternary_tree_ops)
print(cost(new_ternary_tree_ops))

def ternary_tree_test(qubit_op):
    # tests a necessary (but not sufficient) condition for ternary tree
    for i in range(len(qubit_op)):
        for j in range(i, len(qubit_op)):
            iterm = list(qubit_op.terms.items())[i][0]
            jterm = list(qubit_op.terms.items())[j][0]
            count = 0
            for k in range(qubit_op.num_qubits):
                if iterm[k] > 0 and jterm[k] > 0 and iterm[k] != jterm[k]:
                    count += 1
            if count > 1:
                print(iterm, jterm)
                return False
    return True

print(ternary_tree_test(new_ternary_tree_ops))
#print(res[2])
#print(np.linalg.eig(res[2].to_matrix())[0])
