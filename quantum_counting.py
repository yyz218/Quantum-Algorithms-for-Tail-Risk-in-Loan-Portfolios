from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler
import numpy as np
import math
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.circuit.library import Diagonal, GroverOperator
from qiskit.circuit.library import QFT

# quantum counting algorithm to solve the K bits simulation problem

# compute number of solutions M base on the phase
def calculate_M(t, n, hist):
    measured_str = max(hist, key=hist.get)
    measured_int = int(measured_str, 2)
    theta = (measured_int/(2**t))*math.pi*2
    N = 2**n
    M = N * (math.sin(theta/2)**2)
    m = t-1
    err = (2*math.sqrt(M*N) + N/(2**(m+1)))*(2**(-m-1))
    return M, err

# grover's operator
def grover_operator(n, n_iterations, m):
    oracle = Diagonal([-1]*m+[1]*(2**n-m))
    grover_it = GroverOperator(oracle).repeat(n_iterations).to_gate()
    return grover_it


# main function to run the quantum counting algorithm
def main_f(x, n, m):
    qft_dagger = QFT(x, inverse=True).to_gate()
    qft_dagger.label = "QFT†"
    t = x
    qc = QuantumCircuit(n+t, t)

    for qubit in range(t+n):
        qc.h(qubit)

    n_iterations = 1
    for qubit in range(t):
        cgrit = grover_operator(n, n_iterations, m).control()
        qc.append(cgrit, [qubit] + list(range(t, n+t)))
        n_iterations *= 2
    qc.append(qft_dagger, range(t))
    qc.measure(range(t), range(t))
    
    return qc

# simulation
def sim(qc, x, n):
    simulator = AerSimulator()
    transpiled_qc = transpile(qc, simulator)
    job = simulator.run(transpiled_qc, shots=shot)
    hist = job.result().get_counts()
    return calculate_M(x, n, hist)

# actual quantum computer
def actual(qc, service, x, n):
    backend = service.least_busy(operational=True, simulator=False)
    transpiled_qc = transpile(qc, backend=backend)
    job = backend.run([transpiled_qc], shots=shot)
    hist = job.result().get_counts()
    return calculate_M(x, n, hist)

#service = QiskitRuntimeService(channel=channel, instance=instance, token=token)

# set the parameters for the precision x, the total number of elements n
# the number of solutions m, the number of shots for each simulation, and the number of total runs
x, n, m, shot = 6, 5, 6, 10000
qc = main_f(x, n, m)
M, err = sim(qc, x, n)
print(f"The precision = {x}, N = {n}, the actual number of solutions = {m}")
print(f"The number of solutions = {M:.4f}")
print("Error < %.6f" % err)
