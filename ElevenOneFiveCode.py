import qiskit as qs
import numpy as np
import cmath
import math

#"""
output = "1111111111"
stabilizers = [
    [3,3,3,3,3,3,0,0,0,0,0],
    [1,1,1,1,1,1,0,0,0,0,0],
    [0,0,0,3,1,2,2,2,2,1,3],
    [0,0,0,1,2,3,3,3,3,2,1],
    [3,2,1,0,0,0,3,2,1,0,0],
    [1,3,2,0,0,0,1,3,2,0,0],
    [0,0,0,3,2,1,1,2,3,0,0],
    [0,0,0,1,3,2,3,1,2,0,0],
    [3,1,2,0,0,0,3,3,3,1,2],
    [2,3,1,0,0,0,2,2,2,3,1]
]
while(int(output) != 0):
    qc = qs.QuantumCircuit(12, 11)
    for j in range(len(stabilizers)):#does all the stabilizers of the [11,1,5] code
        qc.h(11)
        for i in range(len(stabilizers[j])):
            if(stabilizers[j][i] == 1):
                qc.cx(11,i)
            if(stabilizers[j][i] == 2):
                qc.cy(11,i)
            if(stabilizers[j][i] == 3):
                qc.cz(11,i)
        qc.h(11)
        qc.measure(11,j+1)
        qc.reset(11)
        
    #measures in z basis to get codeword
    qc.h(11)
    qc.cz(11,6)
    qc.cz(11,7)
    qc.cz(11,8)
    qc.cz(11,9)
    qc.cz(11,10)
    qc.h(11)
    qc.measure(11,0)
    qc.reset(11)

    sv_sim = qs.Aer.get_backend('aer_simulator')

    qc.save_statevector()
    sv_job = sv_sim.run(qc,shots=1)
    counts = sv_job.result().get_counts()
    for item in counts:
        output = item
    #print(output)
state = sv_job.result().get_statevector().data

print("Found Zero State")
zeroStatebin = [float(0) for i in range(2**11)]
oneStatebin = []
zeroState = []
oneState = []
for i in range(len(state)):
    if(state[i] > 0):
        #print(i)
        binaryrep = format(i,'b')
        while(len(binaryrep) < 11):
            binaryrep = "0"+binaryrep
        print(binaryrep)
        onebinaryrep = ""
        for j in range(11):
            if(j<6):
                onebinaryrep = onebinaryrep + binaryrep[j]
            else:
                if(binaryrep[j] == "0"):
                    onebinaryrep = onebinaryrep + "1"
                else:
                    onebinaryrep = onebinaryrep + "0"
        zeroStatebin.append(binaryrep)
        zeroState[int(binaryrep, 2)] = 1
        oneStatebin.append(onebinaryrep)
        oneState[int(onebinaryrep, 2)] = 1
        #print(state[i])
print("len(zeroState) = "+str(len(zeroStatebin))+"\n")

print(zeroStatebin)
print(zeroState)
print()
print(oneStatebin)
print(oneState)
#"""
