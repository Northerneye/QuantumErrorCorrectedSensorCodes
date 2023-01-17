#Test the KnillLaflamme condition for error corrected Steane Code sensor
import qiskit as qs
import numpy as np

degree = 0.1

#"""
ZeroLogical = [0 for i in range(2**7)]
ZeroLogical[0] = 1/np.sqrt(8)
ZeroLogical[2**6 + 2**4 + 2**2 + 2**0] = 1/np.sqrt(8)
ZeroLogical[2**5 + 2**4 + 2**1 + 2**0] = 1/np.sqrt(8)
ZeroLogical[2**6 + 2**5 + 2**2 + 2**1] = 1/np.sqrt(8)
ZeroLogical[2**3 + 2**2 + 2**1 + 2**0] = 1/np.sqrt(8)
ZeroLogical[2**6 + 2**4 + 2**3 + 2**1] = 1/np.sqrt(8)
ZeroLogical[2**5 + 2**4 + 2**3 + 2**2] = 1/np.sqrt(8)
ZeroLogical[2**6 + 2**5 + 2**3 + 2**0] = 1/np.sqrt(8)

OneLogical = [0 for i in range(2**7)]
OneLogical[2**7-1] = 1/np.sqrt(8)
OneLogical[2**5 + 2**3 + 2**1] = 1/np.sqrt(8)
OneLogical[2**6 + 2**3 + 2**2] = 1/np.sqrt(8)
OneLogical[2**4 + 2**3 + 2**0] = 1/np.sqrt(8)
OneLogical[2**6 + 2**5 + 2**4] = 1/np.sqrt(8)
OneLogical[2**5 + 2**2 + 2**0] = 1/np.sqrt(8)
OneLogical[2**6 + 2**1 + 2**0] = 1/np.sqrt(8)
OneLogical[2**4 + 2**2 + 2**1] = 1/np.sqrt(8)
"""

ZeroLogical = [0 for i in range(2**7)]
ZeroLogical[0] = 1/np.sqrt(8)
ZeroLogical[2**0 + 2**2 + 2**4 + 2**6] = 1/np.sqrt(8)
ZeroLogical[2**1 + 2**2 + 2**5 + 2**6] = 1/np.sqrt(8)
ZeroLogical[2**0 + 2**1 + 2**4 + 2**5] = 1/np.sqrt(8)
ZeroLogical[2**3 + 2**4 + 2**5 + 2**6] = 1/np.sqrt(8)
ZeroLogical[2**0 + 2**2 + 2**3 + 2**5] = 1/np.sqrt(8)
ZeroLogical[2**1 + 2**2 + 2**3 + 2**4] = 1/np.sqrt(8)
ZeroLogical[2**0 + 2**1 + 2**3 + 2**6] = 1/np.sqrt(8)

OneLogical = [0 for i in range(2**7)]
OneLogical[2**7-1] = 1/np.sqrt(8)
OneLogical[2**1 + 2**3 + 2**5] = 1/np.sqrt(8)
OneLogical[2**0 + 2**3 + 2**4] = 1/np.sqrt(8)
OneLogical[2**2 + 2**3 + 2**6] = 1/np.sqrt(8)
OneLogical[2**0 + 2**1 + 2**2] = 1/np.sqrt(8)
OneLogical[2**1 + 2**4 + 2**6] = 1/np.sqrt(8)
OneLogical[2**0 + 2**5 + 2**6] = 1/np.sqrt(8)
OneLogical[2**2 + 2**4 + 2**5] = 1/np.sqrt(8)


newZeroLogical = ((np.array(ZeroLogical)+np.array(OneLogical))/np.sqrt(2)).tolist()
newOneLogical = ((np.array(ZeroLogical)-np.array(OneLogical))/np.sqrt(2)).tolist()
ZeroLogical = newZeroLogical[:]
OneLogical = newOneLogical[:]
#input(np.matmul(np.matrix(ZeroLogical), np.matrix((OneLogical)).getH()))
#"""

KLZeroStates = []#[ZeroLogical]
KLOneStates = []#[OneLogical]

#########################################################################Getting rotated Z_L|0>
qc = qs.QuantumCircuit(7)
qc.initialize(ZeroLogical, qc.qubits)

CCZPhi = qs.circuit.library.PhaseGate(-4*degree).control(2)
CCZPhiDagger = qs.circuit.library.PhaseGate(4*degree).control(2)

qc.append(CCZPhi, [0,2,4])
qc.append(CCZPhiDagger, [1,3,5])

qc.barrier()
qc.p(degree,0)
qc.p(degree,2)
qc.p(degree,4)
qc.p(degree,6)
qc.barrier()
qc.p(-degree,1)
qc.p(-degree,3)
qc.p(-degree,5)

sv_sim = qs.Aer.get_backend('aer_simulator')

qc.save_statevector()
job = sv_sim.run(qc)

state = job.result().get_statevector().data
#print(state)
#input(qc)
currentState = np.array(state)
KLZeroStates.append(currentState)

##############################################################Getting rotated Z_L|1>
qc = qs.QuantumCircuit(7)
qc.initialize(OneLogical, qc.qubits)

qc.append(CCZPhi, [0,2,4])
qc.append(CCZPhiDagger, [1,3,5])

qc.barrier()
qc.p(degree,0)
qc.p(degree,2)
qc.p(degree,4)
qc.p(degree,6)
qc.barrier()
qc.p(-degree,1)
qc.p(-degree,3)
qc.p(-degree,5)

sv_sim = qs.Aer.get_backend('aer_simulator')

qc.save_statevector()
job = sv_sim.run(qc)

state = job.result().get_statevector().data
#print(state)
#input(qc)
currentState = np.array(state)
KLOneStates.append(currentState)



for i in range(63):
    qc = qs.QuantumCircuit(7)

    qc.initialize(ZeroLogical, qc.qubits)

    if(i<21):
        if(i<7):
            qc.x(i)
        if(i>=7 and i < 14):
            qc.z(i%7)
        if(i>=14):
            qc.y(i%7)


    qc.append(CCZPhi, [0,2,4])
    qc.append(CCZPhiDagger, [1,3,5])

    if(i >= 21 and i < 42):
        if(i<28):
            qc.x(i%7)
        if(i>=28 and i < 35):
            qc.z(i%7)
        if(i>=35):
            qc.y(i%7)

    qc.barrier()
    qc.p(degree,0)
    qc.p(degree,2)
    qc.p(degree,4)
    qc.p(degree,6)
    qc.barrier()
    qc.p(-degree,1)
    qc.p(-degree,3)
    qc.p(-degree,5)

    if(i>=42):
        if(i<49):
            qc.x(i%7)
        if(i>=49 and i < 56):
            qc.z(i%7)
        if(i>=56):
            qc.y(i%7)

    sv_sim = qs.Aer.get_backend('aer_simulator')

    qc.save_statevector()
    job = sv_sim.run(qc)

    state = job.result().get_statevector().data
    #print(state)
    #input(qc)
    currentState = np.array(state)
    KLZeroStates.append(currentState)

for i in range(63):
    qc = qs.QuantumCircuit(7)

    qc.initialize(OneLogical, qc.qubits)

    if(i<21):
        if(i<7):
            qc.x(i)
        if(i>=7 and i < 14):
            qc.z(i%7)
        if(i>=14):
            qc.y(i%7)


    qc.append(CCZPhi, [0,2,4])
    qc.append(CCZPhiDagger, [1,3,5])

    if(i >= 21 and i < 42):
        if(i<28):
            qc.x(i%7)
        if(i>=28 and i < 35):
            qc.z(i%7)
        if(i>=35):
            qc.y(i%7)

    qc.barrier()
    qc.p(degree,0)
    qc.p(degree,2)
    qc.p(degree,4)
    qc.p(degree,6)
    qc.barrier()
    qc.p(-degree,1)
    qc.p(-degree,3)
    qc.p(-degree,5)

    if(i>=42):
        if(i<49):
            qc.x(i%7)
        if(i>=49 and i < 56):
            qc.z(i%7)
        if(i>=56):
            qc.y(i%7)

    sv_sim = qs.Aer.get_backend('aer_simulator')

    qc.save_statevector()
    job = sv_sim.run(qc)

    state = job.result().get_statevector().data
    #print(state)
    #input(qc)
    currentState = np.array(state)
    KLOneStates.append(currentState)

flag = True
l = 0
k = 0
for zerostate in KLZeroStates:
    for onestate in KLOneStates:
        if(np.matmul(np.matrix(zerostate), np.matrix(onestate).getH()) > 0):
            flag = False
            print()
            print("UH OH  l="+str(l)+", k="+str(k))
            input(np.matmul(np.matrix(zerostate), np.matrix(onestate).getH()))
            #print(zerostate)
            #print(onestate)
        k += 1
    l += 1
for l in range(64):
    for m in range(64):
        num1 = np.matmul(np.matrix(KLZeroStates[l]), np.matrix(KLZeroStates[m]).getH()).tolist()[0][0]
        num1 = round(num1.real, 5) + round(num1.imag, 5) * 1j
        num2 = np.matmul(np.matrix(KLOneStates[l]), np.matrix(KLOneStates[m]).getH()).tolist()[0][0]
        num2 = round(num2.real, 5) + round(num2.imag, 5) * 1j
        if(num1 != num2):
            print("Uh Oh!! Unequal "+str(l)+", "+str(m))
            print(str(num1)+" != "+str(num2)+"\n")
#print(len(KLZeroStates))
if(flag):
    print("NONE Found!!")    
input("END")