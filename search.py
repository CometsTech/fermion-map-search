from pauli_list import *
#from queue import PriorityQueue
import heapq
import itertools
import random
from math import log

def search_astar(op, gate_set, timeout=100, debug=0):
    '''
    TODO gate_set should take gates not pairs
    debug: 0 for none, 1 for basic, 2 for more
    '''
    start_cost = cost(op)
    best_cost = start_cost
    best_gates = []
    queue = [(start_cost, [], op)] # priority queue using heapq
    count = 0
    while len(queue) > 0:
        count += 1
        current_cost, gate_keys, current_op = heapq.heappop(queue)
        for next_gate_key in gate_set:
            next_gate = gate_set[next_gate_key]
            next_op = UAUd(next_gate, current_op)
            next_cost = cost(next_op)
            if next_cost < cost(current_op):
                if debug > 1: print(gate_keys + [next_gate], cost(next_op))
                heapq.heappush(queue, (next_cost, gate_keys + [next_gate_key], next_op))
                if next_cost < best_cost:
                    best_cost = next_cost
                    best_gates = gate_keys + [next_gate_key]
                    if debug: print(best_cost, best_gates, count)
                    count = 0
            if debug and np.random.random() < 0.001: print("Random check", count, len(queue))
            if count > timeout:
                break

    best_op = UAUd_list((gate_set[k] for k in best_gates), op)
    return (best_cost, best_gates, best_op)


def search_dfs_CNOT(op, max_depth=5, debug=False):
    start_cost = cost(op)
    best_cost = start_cost
    best_pairs = []
    stack = [[]]
    all_pairs = list(itertools.permutations(range(op.num_qubits), 2))
    while len(stack) > 0:
        pair_list = stack.pop()
        if len(pair_list) > max_depth:
            continue
        current_op = UAUd_CNOT_list(pair_list, op)
        for next_pair in all_pairs:
            next_op = UAUd_CNOT_list([next_pair], current_op)
            next_cost = cost(next_op)
            if next_cost < cost(current_op):
                if debug: print(pair_list + [next_pair], cost(next_op))
                stack.append(pair_list + [next_pair])
                if next_cost < best_cost:
                    best_cost = next_cost
                    best_pairs = pair_list + [next_pair]
    return(best_cost, best_pairs)


def search_SA(op, gate_set, beta_schedule=None, max_depth=3000000, debug=False, store_gates=True):
    start_cost = cost(op)
    current_op = op
    current_cost = cost(current_op)
    
    best_cost=current_cost
    best_op = op
    best_gates=[]

    current_gates = []

    gate_set_keys = list(gate_set.keys())
    beta_schedule = beta_schedule or (lambda n : np.log(1+n/10)/start_cost*50)

    for n in range(max_depth):
        beta = beta_schedule(n)

        next_gate_key = random.choice(gate_set_keys)
        next_op = UAUd(gate_set[next_gate_key], current_op)
        next_cost = cost(next_op)
        prob = np.exp(-beta*(next_cost-current_cost))
        if random.random() < prob:
            if store_gates:
                current_gates = current_gates + [next_gate_key]
            current_cost = next_cost
            current_op = next_op
            
            if current_cost < best_cost:
                best_gates = current_gates
                best_cost = current_cost
                best_op = current_op

            if debug: print(current_cost, best_cost, beta, n, next_gate_key, prob)
            
    return(best_cost, best_gates, current_op)
