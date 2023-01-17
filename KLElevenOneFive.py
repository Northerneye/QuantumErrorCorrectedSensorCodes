import qiskit as qs
import numpy as np
import sympy as sy
import cmath
import math

sy.var('phi x ', real=True)

global myphi
myphi = float(0.1)
#degree = 0.1

myOperator = [[float(0) for i in range(2**11)] for j in range(2**11)]

global Unitary
Unitary = [[expr.subs(x,0) for i in range(2**11)] for j in range(2**11)]


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
zeroStatebin = []
oneStatebin = []
ZeroLogical = [float(0) for i in range(2**11)]
OneLogical = [float(0) for i in range(2**11)]
for i in range(len(state)):
    if(state[i] > 0):
        #print(i)
        binaryrep = format(i,'b')
        while(len(binaryrep) < 11):
            binaryrep = "0"+binaryrep
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
        ZeroLogical[int(binaryrep, 2)] = 1
        oneStatebin.append(onebinaryrep)
        OneLogical[int(onebinaryrep, 2)] = 1
        myOperator[int(binaryrep, 2)][int(binaryrep, 2)] = 1
        Unitary[int(binaryrep, 2)][int(binaryrep, 2)] = 1
        myOperator[int(onebinaryrep, 2)][int(onebinaryrep, 2)] = sy.exp(sy.I*phi)
        Unitary[int(onebinaryrep, 2)][int(onebinaryrep, 2)] = sy.exp(sy.I*phi)
        #print(state[i])

A = sy.IndexedBase('A')
l = 0
for i in range(2**11):
    if(myOperator[i][i] == 0):
        myOperator[i][i] = A[l]
        l += 1
#"""

print("Creating KL conditions...")
#<0|ZxZx|0> - <1|ZxZx|1> = 0
identity = np.array([
    [1,0],
    [0,1]])
X = np.array([
    [0,1],
    [1,0]])
Y = np.array([
    [0, -1*sy.I],
    [sy.I, 0]])
Z = np.array([
    [1, 0],
    [0, -1]])

X0 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(X,identity),identity),identity),identity),identity),identity),identity),identity),identity),identity)
X1 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,X),identity),identity),identity),identity),identity),identity),identity),identity),identity)
X2 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),X),identity),identity),identity),identity),identity),identity),identity),identity)
X3 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),X),identity),identity),identity),identity),identity),identity),identity)
X4 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),X),identity),identity),identity),identity),identity),identity)
X5 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),X),identity),identity),identity),identity),identity)
X6 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),X),identity),identity),identity),identity)
X7 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),identity),X),identity),identity),identity)
X8 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),identity),identity),X),identity),identity)
X9 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),identity),identity),identity),X),identity)
X10 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),identity),identity),identity),identity),X)
Y0 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(Y,identity),identity),identity),identity),identity),identity),identity),identity),identity),identity)
Y1 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,Y),identity),identity),identity),identity),identity),identity),identity),identity),identity)
Y2 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),Y),identity),identity),identity),identity),identity),identity),identity),identity)
Y3 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),Y),identity),identity),identity),identity),identity),identity),identity)
Y4 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),Y),identity),identity),identity),identity),identity),identity)
Y5 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),Y),identity),identity),identity),identity),identity)
Y6 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),Y),identity),identity),identity),identity)
Y7 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),identity),Y),identity),identity),identity)
Y8 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),identity),identity),Y),identity),identity)
Y9 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),identity),identity),identity),Y),identity)
Y10 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),identity),identity),identity),identity),Y)
Z0 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(Z,identity),identity),identity),identity),identity),identity),identity),identity),identity),identity)
Z1 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,Z),identity),identity),identity),identity),identity),identity),identity),identity),identity)
Z2 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),Z),identity),identity),identity),identity),identity),identity),identity),identity)
Z3 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),Z),identity),identity),identity),identity),identity),identity),identity)
Z4 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),Z),identity),identity),identity),identity),identity),identity)
Z5 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),Z),identity),identity),identity),identity),identity)
Z6 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),Z),identity),identity),identity),identity)
Z7 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),identity),Z),identity),identity),identity)
Z8 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),identity),identity),Z),identity),identity)
Z9 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),identity),identity),identity),Z),identity)
Z10 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),identity),identity),identity),identity),Z)



#complete set of KL conditions (still only find 103 conditions total)
expr = x
Ab = []
Elist = [X0]#, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9, Y10, Z0, Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9, Z10]
Elistname = ["X0"]#, "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "Y0", "Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9", "Y10", "Z0", "Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z8", "Z9", "Z10"]
localcounter = 0
for error1 in Elist:
    for error2 in Elist:
        errorinnerproduct = np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(myOperator),np.matrix(error1))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error2),np.matrix(myOperator))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(myOperator),np.matrix(error1))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error2),np.matrix(myOperator))).getH()).tolist()[0][0]
        Ab.append(errorinnerproduct)
        print(errorinnerproduct)
        #input(Ab)
        errorinnerproduct = np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error1),np.matrix(myOperator))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error2),np.matrix(myOperator))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error1),np.matrix(myOperator))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error2),np.matrix(myOperator))).getH()).tolist()[0][0]
        Ab.append(errorinnerproduct)
        print(errorinnerproduct)
        #input(Ab)
        errorinnerproduct = np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(myOperator),np.matrix(error1))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(myOperator),np.matrix(error2))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(myOperator),np.matrix(error1))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(myOperator),np.matrix(error2))).getH()).tolist()[0][0]
        Ab.append(errorinnerproduct)
        print(errorinnerproduct)
        #input(Ab)
        print("Completed Error "+str(Elistname[localcounter%len(Elist)])+" "+str(Elistname[math.floor(localcounter/len(Elist))]))
        localcounter += 1
#"""

#Ab = sy.Tuple(i for i in Ab)
Ab = sy.Tuple.fromiter(Ab)

#SYMPY SOLVER MAY TAKE A LONG TIME TO FINISH
solution = sy.solve(Ab,Ab.atoms(sy.Indexed))
print(solution)

global currentmatrix

theoretical = True
l = 0
for i in range(2**11):
    if(Unitary[i][i] == 0):
        Unitary[i][i] = solution.args[0][l]#.subs([(var101, -var108),(var102, -var107),(var103,-var106)])#SEE EXPRESSION FOR GATES as a function of free variables CCZ(var110,var109...)
        l += 1
Unitary = sy.Matrix(Unitary)
identity = np.array([
    [1,0],
    [0,1]])
currentmatrix = sy.Matrix(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),identity))

global chargates 
chargates = []

def makeqiskitgate(curbinary, i):
    global currentmatrix
    global Unitary
    global chargates
    global myphi
    global Unitarysubbed
    mybits = []
    identity = np.array([
        [1,0],
        [0,1]])
    while(len(curbinary) < 11):
        curbinary = "0"+curbinary
    zflag = False
    workingmatrix = np.array([1])
    for j in range(len(curbinary)):
        if(curbinary[j] == "1"):
            if(zflag==False):
                #need to make current diagonal entry equal to desired
                nextmatrix = np.array([[0, 0], #Z((Unitary*currentmatrix.H)[i,i])
                 [0, 1-sy.simplify((Unitarysubbed*currentmatrix.H)[i,i])]])
                workingmatrix = np.kron(workingmatrix,nextmatrix)
                zflag = True
                chargates.append("(Z_"+str(j)+"="+str(sy.simplify((Unitarysubbed*currentmatrix.H)[i,i]))+")")
                mybits.append(j)
            else: #kronecker product a control signal
                nextmatrix = np.array([
                    [0,0],
                    [0,1]])
                workingmatrix = np.kron(workingmatrix, nextmatrix)
                chargates[len(chargates)-1] = chargates[len(chargates)-1] + "C_"+str(j)
                mybits.append(j)
        else: # kronecker product an identity
            workingmatrix = np.kron(workingmatrix, identity)
    
    workingmatrix = sy.Matrix(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),identity),identity),identity),identity),identity)) - sy.Matrix(workingmatrix)
    
    if(len(mybits) > 1):
        nextgate = qs.circuit.library.PhaseGate(float(sy.simplify(-sy.I*sy.ln((Unitarysubbed*currentmatrix.H)[i,i])).subs(phi,myphi).evalf())).control(len(mybits)-1)
    else:
        nextgate = qs.circuit.library.PhaseGate(float(sy.simplify(-sy.I*sy.ln((Unitarysubbed*currentmatrix.H)[i,i])).subs(phi,myphi).evalf()))

    #return identity minus calculated matrix   
    #return sy.Matrix(np.kron(identity,identity)) - sy.Matrix(workingmatrix)
    
    # 7 qubit return   
    return [nextgate, mybits, sy.Matrix(workingmatrix)]

def makegate(curbinary, i):
    global currentmatrix
    global Unitary
    global chargates
    global myphi
    identity = np.array([
        [1,0],
        [0,1]])
    while(len(curbinary) < 11):
        curbinary = "0"+curbinary
    zflag = False
    workingmatrix = np.array([1])
    for j in range(len(curbinary)):
        if(curbinary[j] == "1"):
            if(zflag==False):
                #need to make current diagonal entry equal to desired
                nextmatrix = np.array([[0, 0], #Z((Unitary*currentmatrix.H)[i,i])
                 [0, 1-sy.simplify((Unitary*currentmatrix.H)[i,i])]])
                workingmatrix = np.kron(workingmatrix,nextmatrix)
                zflag = True
                chargates.append("(Z_"+str(j)+"="+str(sy.simplify((Unitary*currentmatrix.H)[i,i]))+")")
            else: #kronecker product a control signal
                nextmatrix = np.array([
                    [0,0],
                    [0,1]])
                workingmatrix = np.kron(workingmatrix, nextmatrix)
                chargates[len(chargates)-1] = chargates[len(chargates)-1] + "C_"+str(j)
        else: # kronecker product an identity
            workingmatrix = np.kron(workingmatrix, identity)
    

    
    #return identity minus calculated matrix   
    #return sy.Matrix(np.kron(identity,identity)) - sy.Matrix(workingmatrix)
    
    # 7 qubit return   
    return sy.Matrix(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),identity),identity),identity),identity),identity)) - sy.Matrix(workingmatrix)
    


qc = qs.QuantumCircuit(11)
print("ASSEMBLING CIRCUIT...")
global Unitarysubbed
#NEED TO SUB WITH INDEXED VARIABLES INSTEAD (A[1] instead of var10)
Unitarysubbed = Unitary.subs([(var101, -var108),(var102, -var107),(var103,-var106)]).subs([(var98, -1),(var99, 1),(var100, -1),(var106, -1),(var107, 1),(var108, -1),(var109, 1),(var110, -1),(var111, 1)]).evalf()
for i in range(2**11):
    if(Unitarysubbed[i,i] != currentmatrix[i,i]):
        #Use i to find necessary gate(from binary decomp)(01101 is c2c3z5)
        curbinary = format(i,'b')
        [nextgate, mybits, newMatrix] = makeqiskitgate(curbinary, i)
        qc.append(nextgate, mybits)
        currentmatrix = currentmatrix * newMatrix
        
print(chargates)
print("gates(phi) ^^\n")



#Elist = [X0, X1, X2, X3, X4, X5, X6, Y0, Y1, Y2, Y3, Y4, Y5, Y6, Z0, Z1, Z2, Z3, Z4, Z5, Z6]
#Elistname = ["X0", "X1", "X2", "X3", "X4", "X5", "X6", "Y0", "Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Z0", "Z1", "Z2", "Z3", "Z4", "Z5", "Z6"]
Elist = [X3]
Elistname = ["X3"]
localcounter = 0
for error1 in Elist:
    for error2 in Elist:
        print()
        print(np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(myOperator),np.matrix(error1))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error2),np.matrix(myOperator))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(myOperator),np.matrix(error1))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error2),np.matrix(myOperator))).getH()).tolist()[0][0])
        print(np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(Unitary),np.matrix(error1))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error2),np.matrix(Unitary))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(Unitary),np.matrix(error1))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error2),np.matrix(Unitary))).getH()).tolist()[0][0])
        print(np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(Unitarysubbed),np.matrix(error1))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error2),np.matrix(Unitarysubbed))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(Unitarysubbed),np.matrix(error1))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error2),np.matrix(Unitarysubbed))).getH()).tolist()[0][0])
        print()
        print(np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error1),np.matrix(myOperator))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error2),np.matrix(myOperator))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error1),np.matrix(myOperator))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error2),np.matrix(myOperator))).getH()).tolist()[0][0])
        print(np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error1),np.matrix(Unitary))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error2),np.matrix(Unitary))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error1),np.matrix(Unitary))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error2),np.matrix(Unitary))).getH()).tolist()[0][0])
        print(np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error1),np.matrix(Unitarysubbed))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error2),np.matrix(Unitarysubbed))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error1),np.matrix(Unitarysubbed))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error2),np.matrix(Unitarysubbed))).getH()).tolist()[0][0])
        print()
        print(np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(myOperator),np.matrix(error1))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(myOperator),np.matrix(error2))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(myOperator),np.matrix(error1))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(myOperator),np.matrix(error2))).getH()).tolist()[0][0])
        print(np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(Unitary),np.matrix(error1))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(Unitary),np.matrix(error2))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(Unitary),np.matrix(error1))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(Unitary),np.matrix(error2))).getH()).tolist()[0][0])
        print(np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(Unitarysubbed),np.matrix(error1))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(Unitarysubbed),np.matrix(error2))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(Unitarysubbed),np.matrix(error1))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(Unitarysubbed),np.matrix(error2))).getH()).tolist()[0][0])
        print("Error X3X3 not "+str(Elistname[localcounter%len(Elist)])+" "+str(Elistname[math.floor(localcounter/len(Elist))])+"\n")
        localcounter += 1
input("ALL KL Condtions Calculated")

print(qc)
print("qiskit gates with phi = "+str(myphi))

mygate = qs.quantum_info.Operator(np.array(Unitarysubbed.subs(phi,myphi)))


qubitorder = [0,1,2,3,4,5,6,7,8,9,10]

KLZeroStates = []#[ZeroLogical]
KLOneStates = []#[OneLogical]

logicalStates = [ZeroLogical, OneLogical]
for j in range(2):
    for i in range(66):
        qc = qs.QuantumCircuit(11)

        qc.initialize(logicalStates[j], qc.qubits)

        if(i<33):
            if(i<11):
                qc.x(i)
            if(i>=11 and i < 22):
                qc.z(i%11)
            if(i>=33):
                qc.y(i%11)

        qc.unitary(mygate, qubitorder)

        if(i >= 44 and i < 66):
            if(i<44):
                qc.x(i%11)
            if(i>=44 and i < 55):
                qc.z(i%11)
            if(i>=66):
                qc.y(i%11)

        sv_sim = qs.Aer.get_backend('aer_simulator')

        qc.save_statevector()
        job = sv_sim.run(qc)

        state = job.result().get_statevector().data
        #print(state)
        #input(qc)
        currentState = np.array(state)
        if(j==0):
            KLZeroStates.append(currentState)
        elif(j==1):
            KLOneStates.append(currentState)
    #########################################################################Getting rotated Z_L|0>
    qc = qs.QuantumCircuit(11)
    qc.initialize(logicalStates[j], qc.qubits)

    qc.unitary(mygate, qubitorder)

    sv_sim = qs.Aer.get_backend('aer_simulator')

    qc.save_statevector()
    job = sv_sim.run(qc)

    state = job.result().get_statevector().data
    #print(state)
    #input(qc)
    currentState = np.array(state)
    if(j==0):
        KLZeroStates.append(currentState)
    elif(j==1):
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
for l in range(len(KLZeroStates)):
    for m in range(len(KLOneStates)):
        num1 = np.matmul(np.matrix(KLZeroStates[l]), np.matrix(KLZeroStates[m]).getH()).tolist()[0][0]
        num1 = round(num1.real, 5) + round(num1.imag, 5) * 1j
        num2 = np.matmul(np.matrix(KLOneStates[l]), np.matrix(KLOneStates[m]).getH()).tolist()[0][0]
        num2 = round(num2.real, 5) + round(num2.imag, 5) * 1j
        if(num1 != num2):
            print("Uh Oh!! Unequal "+str(l)+", "+str(m))
            print(str(num1)+" != "+str(num2)+"\n")
            flag = False
#print(len(KLZeroStates))
if(flag):
    print("NONE Found!!")  
else:
    print("FAILED")  
print("END")