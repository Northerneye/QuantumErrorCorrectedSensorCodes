import qiskit as qs
import numpy as np
import sympy as sy
import cmath
import math

myvars = ''
for i in range(112):
    myvars = myvars + ' var'+str(i)
sy.var('phi x '+ myvars, real=True)
myvar = [var0,var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11,var12,var13,var14,var15,var16,var17,var18,var19,var20,var21,var22,var23,var24,var25,var26,var27,var28,var29,var30,var31,var32,var33,var34,var35,var36,var37,var38,var39,var40,var41,var42,var43,var44,var45,var46,var47,var48,var49,var50,var51,var52,var53,var54,var55,var56,var57,var58,var59,var60,var61,var62,var63,var64,var65,var66,var67,var68,var69,var70,var71,var72,var73,var74,var75,var76,var77,var78,var79,var80,var81,var82,var83,var84,var85,var86,var87,var88,var89,var90,var91,var92,var93,var94,var95,var96,var97,var98,var99,var100,var101,var102,var103,var104,var105,var106,var107,var108,var109,var110,var111]

global myphi
myphi = float(0.1)
#degree = 0.1

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



myOperator = [[float(0) for i in range(2**7)] for j in range(2**7)]

myOperator[0][0] = 1
myOperator[2**6 + 2**4 + 2**2 + 2**0][2**6 + 2**4 + 2**2 + 2**0] = 1
myOperator[2**5 + 2**4 + 2**1 + 2**0][2**5 + 2**4 + 2**1 + 2**0] = 1
myOperator[2**6 + 2**5 + 2**2 + 2**1][2**6 + 2**5 + 2**2 + 2**1] = 1
myOperator[2**3 + 2**2 + 2**1 + 2**0][2**3 + 2**2 + 2**1 + 2**0] = 1
myOperator[2**6 + 2**4 + 2**3 + 2**1][2**6 + 2**4 + 2**3 + 2**1] = 1
myOperator[2**5 + 2**4 + 2**3 + 2**2][2**5 + 2**4 + 2**3 + 2**2] = 1
myOperator[2**6 + 2**5 + 2**3 + 2**0][2**6 + 2**5 + 2**3 + 2**0] = 1

myOperator[2**7-1][2**7-1] = sy.exp(sy.I*phi)
myOperator[2**5 + 2**3 + 2**1][2**5 + 2**3 + 2**1] = sy.exp(sy.I*phi)
myOperator[2**6 + 2**3 + 2**2][2**6 + 2**3 + 2**2] = sy.exp(sy.I*phi)
myOperator[2**4 + 2**3 + 2**0][2**4 + 2**3 + 2**0] = sy.exp(sy.I*phi)
myOperator[2**6 + 2**5 + 2**4][2**6 + 2**5 + 2**4] = sy.exp(sy.I*phi)
myOperator[2**5 + 2**2 + 2**0][2**5 + 2**2 + 2**0] = sy.exp(sy.I*phi)
myOperator[2**6 + 2**1 + 2**0][2**6 + 2**1 + 2**0] = sy.exp(sy.I*phi)
myOperator[2**4 + 2**2 + 2**1][2**4 + 2**2 + 2**1] = sy.exp(sy.I*phi)

l = 0
for i in range(2**7):
    if(myOperator[i][i] == 0):
        myOperator[i][i] = myvar[l]
        l += 1
myOperator = sy.Matrix(myOperator)

#TEST IF THE GATE COMMUTES WITH STABILIZERS
stabilizersData = np.array([
    [0,0,0,3,3,3,3],
    [0,3,3,0,0,3,3],
    [3,0,3,0,3,0,3],
    [0,0,0,1,1,1,1],
    [0,1,1,0,0,1,1],
    [1,0,1,0,1,0,1]
])
ZeroMatrix = np.array([[float(0) for i in range(2**7)] for j in range(2**7)])
stabilizers = []
for i in range(6):
    qcStab = qs.QuantumCircuit(7)
    for j in range(len(stabilizersData[i])):
        if(stabilizersData[i][j] == 1):
            qcStab.x(j)
        if(stabilizersData[i][j] == 2):
            qcStab.y(j)
        if(stabilizersData[i][j] == 3):
            qcStab.z(j)
    #input(qcStab)
    #input(qs.quantum_info.operators.Operator(qcStab.to_gate()).data)
    mystabilizer = qs.quantum_info.operators.Operator(qcStab.to_gate()).data.tolist()
    for j in range(len(mystabilizer)):
        for k in range(len(mystabilizer[j])):
            mystabilizer[j][k] = mystabilizer[j][k].real
    stabilizers.append(sy.Matrix(mystabilizer))
    
expr = x
A = sy.Matrix([[float(0) for i in range(112)] for j in range(352)])
b = sy.Matrix([float(0) for j in range(352)])
Ab = [[expr.subs(x,0) for i in range(113)] for j in range(1)]
l = 0
o = 0
for stabilizer in stabilizers:
    l += 1
    localcounter = 0
    commutator = myOperator*stabilizer - stabilizer*myOperator
    for i in range(2**7):
        for j in range(2**7):
            if(commutator[i,j] != 0):
                localcounter += 1
                Ab.append([expr.subs(x,0) for i in range(113)])
                for k in range(112):#find coeff of all variables
                    Ab[o][k] = commutator[i,j].coeff(myvar[k])
                Ab[o][112] = commutator[i,j].subs([(var0, 0),(var1, 0),(var2, 0),(var3, 0),(var4, 0),(var5, 0),(var6, 0),(var7, 0),(var8, 0),(var9, 0),(var10, 0),(var11, 0),(var12, 0),(var13, 0),(var14, 0),(var15, 0),(var16, 0),(var17, 0),(var18, 0),(var19, 0),(var20, 0),(var21, 0),(var22, 0),(var23, 0),(var24, 0),(var25, 0),(var26, 0),(var27, 0),(var28, 0),(var29, 0),(var30, 0),(var31, 0),(var32, 0),(var33, 0),(var34, 0),(var35, 0),(var36, 0),(var37, 0),(var38, 0),(var39, 0),(var40, 0),(var41, 0),(var42, 0),(var43, 0),(var44, 0),(var45, 0),(var46, 0),(var47, 0),(var48, 0),(var49, 0),(var50, 0),(var51, 0),(var52, 0),(var53, 0),(var54, 0),(var55, 0),(var56, 0),(var57, 0),(var58, 0),(var59, 0),(var60, 0),(var61, 0),(var62, 0),(var63, 0),(var64, 0),(var65, 0),(var66, 0),(var67, 0),(var68, 0),(var69, 0),(var70, 0),(var71, 0),(var72, 0),(var73, 0),(var74, 0),(var75, 0),(var76, 0),(var77, 0),(var78, 0),(var79, 0),(var80, 0),(var81, 0),(var82, 0),(var83, 0),(var84, 0),(var85, 0),(var86, 0),(var87, 0),(var88, 0),(var89, 0),(var90, 0),(var91, 0),(var92, 0),(var93, 0),(var94, 0),(var95, 0),(var96, 0),(var97, 0),(var98, 0),(var99, 0),(var100, 0),(var101, 0),(var102, 0),(var103, 0),(var104, 0),(var105, 0),(var106, 0),(var107, 0),(var108, 0),(var109, 0),(var110, 0),(var111, 0)]).evalf()
                o += 1
                #input("problem, "+str(commutator[i,j])+" != 0,     i="+str(i)+", j="+str(j)+"\n")
    print("finished stabilizer "+str(l))
    print("o = "+str(o))
    if(localcounter):
        print("Logical Rotation Does Not Commute With Operator "+str(l)+"\n"+str(stabilizersData[l-1]))
        print(str(localcounter)+" equations\n")

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

X0 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(X,identity),identity),identity),identity),identity),identity)
X1 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,X),identity),identity),identity),identity),identity)
X2 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),X),identity),identity),identity),identity)
X3 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),X),identity),identity),identity)
X4 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),X),identity),identity)
X5 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),X),identity)
X6 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),X)
Y0 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(Y,identity),identity),identity),identity),identity),identity)
Y1 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,Y),identity),identity),identity),identity),identity)
Y2 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),Y),identity),identity),identity),identity)
Y3 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),Y),identity),identity),identity)
Y4 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),Y),identity),identity)
Y5 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),Y),identity)
Y6 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),Y)
Z0 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(Z,identity),identity),identity),identity),identity),identity)
Z1 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,Z),identity),identity),identity),identity),identity)
Z2 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),Z),identity),identity),identity),identity)
Z3 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),Z),identity),identity),identity)
Z4 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),Z),identity),identity)
Z5 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),Z),identity)
Z6 = np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),Z)


Elist = [X0, X1, X2, X3, X4, X5, X6]#, Y0, Y1, Y2, Y3, Y4, Y5, Y6, Z0, Z1, Z2, Z3, Z4, Z5, Z6]
Elistname = ["X0", "X1", "X2", "X3", "X4", "X5", "X6", "Y0", "Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Z0", "Z1", "Z2", "Z3", "Z4", "Z5", "Z6"]
if(False):#makes KL conditions into linear system of equations to be solved later(likely nonlinear, ex. var101**2 will cause a failure)
    localcounter = 0
    for error in Elist:
        errorinnerproduct = np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(myOperator),np.matrix(error))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error),np.matrix(myOperator))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(myOperator),np.matrix(error))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error),np.matrix(myOperator))).getH()).tolist()[0][0]
        Ab.append([expr.subs(x,0) for i in range(113)])
        for k in range(112):#find coeff of all variables
            Ab[o][k] = errorinnerproduct.coeff(myvar[k])             
        Ab[o][112] = 0
        o += 1
        errorinnerproduct = np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error),np.matrix(myOperator))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error),np.matrix(myOperator))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error),np.matrix(myOperator))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error),np.matrix(myOperator))).getH()).tolist()[0][0]
        Ab.append([expr.subs(x,0) for i in range(113)])
        for k in range(112):#find coeff of all variables
            Ab[o][k] = errorinnerproduct.coeff(myvar[k])             
        Ab[o][112] = 0
        o += 1
        errorinnerproduct = np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(myOperator),np.matrix(error))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(myOperator),np.matrix(error))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(myOperator),np.matrix(error))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(myOperator),np.matrix(error))).getH()).tolist()[0][0]
        Ab.append([expr.subs(x,0) for i in range(113)])
        for k in range(112):#find coeff of all variables
            Ab[o][k] = errorinnerproduct.coeff(myvar[k])             
        Ab[o][112] = 0
        print("Completed Error "+str(Elistname[localcounter%len(Elist)]))
        localcounter += 1
        o += 1


"""#complete set of KL conditions (still only find 103 conditions total)
Elist = [X0, X1, X2, X3, X4, X5, X6, Y0, Y1, Y2, Y3, Y4, Y5, Y6, Z0, Z1, Z2, Z3, Z4, Z5, Z6]
Elistname = ["X0", "X1", "X2", "X3", "X4", "X5", "X6", "Y0", "Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Z0", "Z1", "Z2", "Z3", "Z4", "Z5", "Z6"]
localcounter = 0
for error1 in Elist:
    for error2 in Elist:
        errorinnerproduct = np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(myOperator),np.matrix(error1))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error2),np.matrix(myOperator))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(myOperator),np.matrix(error1))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error2),np.matrix(myOperator))).getH()).tolist()[0][0]
        Ab.append([expr.subs(x,0) for i in range(113)])
        for k in range(112):#find coeff of all variables
            Ab[o][k] = errorinnerproduct.coeff(myvar[k])             
        Ab[o][112] = 0
        o += 1
        errorinnerproduct = np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error1),np.matrix(myOperator))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(error2),np.matrix(myOperator))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error1),np.matrix(myOperator))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(error2),np.matrix(myOperator))).getH()).tolist()[0][0]
        Ab.append([expr.subs(x,0) for i in range(113)])
        for k in range(112):#find coeff of all variables
            Ab[o][k] = errorinnerproduct.coeff(myvar[k])             
        Ab[o][112] = 0
        o += 1
        errorinnerproduct = np.matmul(np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(myOperator),np.matrix(error1))),np.matmul(np.matrix(ZeroLogical),np.matmul(np.matrix(myOperator),np.matrix(error2))).getH()).tolist()[0][0] - np.matmul(np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(myOperator),np.matrix(error1))),np.matmul(np.matrix(OneLogical),np.matmul(np.matrix(myOperator),np.matrix(error2))).getH()).tolist()[0][0]
        Ab.append([expr.subs(x,0) for i in range(113)])
        for k in range(112):#find coeff of all variables
            Ab[o][k] = errorinnerproduct.coeff(myvar[k])             
        Ab[o][112] = 0
        print("Completed Error "+str(Elistname[localcounter%len(Elist)])+" "+str(Elistname[math.floor(localcounter/len(Elist))]))
        localcounter += 1
        o += 1
#"""

_, inds = sy.Matrix(Ab).T.rref()
print(str(len(inds))+" independent equations")
solution = sy.linsolve(sy.Matrix(np.array(Ab)[list(inds)]), (var0,var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11,var12,var13,var14,var15,var16,var17,var18,var19,var20,var21,var22,var23,var24,var25,var26,var27,var28,var29,var30,var31,var32,var33,var34,var35,var36,var37,var38,var39,var40,var41,var42,var43,var44,var45,var46,var47,var48,var49,var50,var51,var52,var53,var54,var55,var56,var57,var58,var59,var60,var61,var62,var63,var64,var65,var66,var67,var68,var69,var70,var71,var72,var73,var74,var75,var76,var77,var78,var79,var80,var81,var82,var83,var84,var85,var86,var87,var88,var89,var90,var91,var92,var93,var94,var95,var96,var97,var98,var99,var100,var101,var102,var103,var104,var105,var106,var107,var108,var109,var110,var111))
print(solution)

global currentmatrix
global Unitary

Unitary = [[expr.subs(x,0) for i in range(2**7)] for j in range(2**7)]

Unitary[0][0] = 1
Unitary[2**6 + 2**4 + 2**2 + 2**0][2**6 + 2**4 + 2**2 + 2**0] = 1
Unitary[2**5 + 2**4 + 2**1 + 2**0][2**5 + 2**4 + 2**1 + 2**0] = 1
Unitary[2**6 + 2**5 + 2**2 + 2**1][2**6 + 2**5 + 2**2 + 2**1] = 1
Unitary[2**3 + 2**2 + 2**1 + 2**0][2**3 + 2**2 + 2**1 + 2**0] = 1
Unitary[2**6 + 2**4 + 2**3 + 2**1][2**6 + 2**4 + 2**3 + 2**1] = 1
Unitary[2**5 + 2**4 + 2**3 + 2**2][2**5 + 2**4 + 2**3 + 2**2] = 1
Unitary[2**6 + 2**5 + 2**3 + 2**0][2**6 + 2**5 + 2**3 + 2**0] = 1

Unitary[2**7-1][2**7-1] = sy.exp(sy.I*phi)
Unitary[2**5 + 2**3 + 2**1][2**5 + 2**3 + 2**1] = sy.exp(sy.I*phi)
Unitary[2**6 + 2**3 + 2**2][2**6 + 2**3 + 2**2] = sy.exp(sy.I*phi)
Unitary[2**4 + 2**3 + 2**0][2**4 + 2**3 + 2**0] = sy.exp(sy.I*phi)
Unitary[2**6 + 2**5 + 2**4][2**6 + 2**5 + 2**4] = sy.exp(sy.I*phi)
Unitary[2**5 + 2**2 + 2**0][2**5 + 2**2 + 2**0] = sy.exp(sy.I*phi)
Unitary[2**6 + 2**1 + 2**0][2**6 + 2**1 + 2**0] = sy.exp(sy.I*phi)
Unitary[2**4 + 2**2 + 2**1][2**4 + 2**2 + 2**1] = sy.exp(sy.I*phi)

theoretical = True
l = 0
for i in range(2**7):
    if(Unitary[i][i] == 0):
        if(theoretical):#examine solution(the variable) to see, varX**2 ==1 => subs(...)
            Unitary[i][i] = solution.args[0][l]#.subs([(var101, -var108),(var102, -var107),(var103,-var106)])#SEE EXPRESSION FOR GATES as a function of free variables CCZ(var110,var109...)
        else:#See gates for a specific value of the free variables
            Unitary[i][i] = solution.args[0][l].subs([(var101, -var108),(var102, -var107),(var103,-var106)]).subs([(var98, 1),(var99, 1),(var100, 1),(var106, 1),(var107, 1),(var108, 1),(var109, 1),(var110, 1),(var111, 1)]).evalf()
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
    while(len(curbinary) < 7):
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
    
    workingmatrix = sy.Matrix(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),identity)) - sy.Matrix(workingmatrix)
    
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
    while(len(curbinary) < 7):
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
    return sy.Matrix(np.kron(np.kron(np.kron(np.kron(np.kron(np.kron(identity,identity),identity),identity),identity),identity),identity)) - sy.Matrix(workingmatrix)


"""#BRUTE FORCE CHECKS THE FREE VARIABLES
myphases = [1, phi]
qc = qs.QuantumCircuit(7)
varleft = [var98, var99, var100, var101, var102, var103, var106, var107, var108, var109, var110, var111]
numgates = 1000
numgatesindex = [0 for i in range(len(varleft))]
print("Searching for smallest circuit...")
for j in range(len(myphases)**len(varleft)):
    mybinary = format(j, 'b')
    while(len(mybinary) < len(varleft)):
        mybinary = "0"+mybinary
    localUnitary = Unitary.subs([(varleft[0], myphases[int(mybinary[0])]),(varleft[1], myphases[int(mybinary[1])]),(varleft[2], myphases[int(mybinary[2])]),(varleft[3], myphases[int(mybinary[3])]),(varleft[4], myphases[int(mybinary[4])]),(varleft[5], myphases[int(mybinary[5])]),(varleft[6], myphases[int(mybinary[6])]),(varleft[7], myphases[int(mybinary[7])]),(varleft[8], myphases[int(mybinary[8])]),(varleft[9], myphases[int(mybinary[9])]),(varleft[10], myphases[int(mybinary[10])]),(varleft[11], myphases[int(mybinary[11])]),])
    print(mybinary)
    localnumgates = 0
    for i in range(2**7):
        if(Unitary[i,i] != currentmatrix[i,i]):
            #Use i to find necessary gate(from binary decomp)(01101 is c2c3z5)
            curbinary = format(i,'b')
            newMatrix = makegate(curbinary, i)
            currentmatrix = currentmatrix * newMatrix
            localnumgates += 1
    if(localnumgates < numgates):
        print(localnumgates)
        numgates = localnumgates
        numgatesindex = mybinary
            
print(curbinary)
#"""

qc = qs.QuantumCircuit(7)
print("ASSEMBLING CIRCUIT...")
global Unitarysubbed
Unitarysubbed = Unitary.subs([(var101, -var108),(var102, -var107),(var103,-var106)]).subs([(var98, -1),(var99, 1),(var100, -1),(var106, -1),(var107, 1),(var108, -1),(var109, 1),(var110, -1),(var111, 1)]).evalf()
for i in range(2**7):
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

"""#example of gate decomposition algorithm
currentmatrix = sy.Matrix(np.kron(identity,identity))

Unitary = sy.Matrix([
    [1,  0,  0,  0],
    [0,  sy.exp(sy.I*a),  0,  0],
    [0,  0,  sy.exp(sy.I*a),  0],
    [0,  0,  0, sy.exp(sy.I*a)]
])

for i in range(2**2):
    if(Unitary[i,i] != currentmatrix[i,i]):
        #Use i to find necessary gate(from binary decomp)(01101 is c2c3z5)
        curbinary = format(i,'b')
        newMatrix = makegate(curbinary, i)
        print()
        print(newMatrix)
        currentmatrix = currentmatrix * newMatrix
        print(currentmatrix)
        print(chargates)
        input(Unitary)
input("finished")
"""


"""
l = 0
for stabilizer in stabilizers:
    l += 1
    localflag = False
    commutator = np.matmul(np.matrix(mygate), np.matrix(stabilizer)) - np.matmul(np.matrix(stabilizer), np.matrix(mygate)).tolist()
    for i in range(len(commutator)):
        for j in range(len(commutator[i].tolist()[0])):
            if(round(commutator[i].tolist()[0][j].real, 5) != 0 and round(commutator[i].tolist()[0][j].imag, 5) != 0):
                localflag = True
                #print("problem, "+str(commutator[i].tolist()[0][j])+" != 0,     i="+str(i)+", j="+str(j)+"\n")
    if(localflag):
        print("Logical Rotation Does Not Commute With Operator "+str(l)+"\n"+str(stabilizersData[l-1]))
#"""
"""#TEST OVERLAPS AND OUTPUT STATE
qc.save_statevector()
sv_job = sv_sim.run(qc)
state = sv_job.result().get_statevector().data
ZeroOverlap = np.matmul(np.matrix(np.array(ZeroLogical)), np.matrix(np.array(state)).getH()) #<0|psi>
OneOverlap = np.matmul(np.matrix(np.array(OneLogical)), np.matrix(np.array(state)).getH()) #<1|psi>
PlusOverlap = np.matmul(np.matrix((np.array(ZeroLogical)+np.array(OneLogical))/np.sqrt(2)), np.matrix(np.array(state)).getH()) #<+|psi>
MinusOverlap = np.matmul(np.matrix((np.array(ZeroLogical)-np.array(OneLogical))/np.sqrt(2)), np.matrix(np.array(state)).getH()) #<-|psi>
PlusIOverlap = np.matmul(np.matrix((np.array(ZeroLogical)+np.array(OneLogical)*1j)/np.sqrt(2)), np.matrix(np.array(state)).getH()) #<+i|psi>
MinusIOverlap = np.matmul(np.matrix((np.array(ZeroLogical)-np.array(OneLogical)*1j)/np.sqrt(2)), np.matrix(np.array(state)).getH()) #<-i|psi>
print("|<0|psi>|^2 = "+str(ZeroOverlap*np.conjugate(ZeroOverlap))+"\n")
print("|<1|psi>|^2 = "+str(OneOverlap*np.conjugate(OneOverlap))+"\n")
print("|<+|psi>|^2 = "+str(PlusOverlap*np.conjugate(PlusOverlap))+"\n")
print("|<-|psi>|^2 = "+str(MinusOverlap*np.conjugate(MinusOverlap))+"\n")
print("|<+i|psi>|^2 = "+str(PlusIOverlap*np.conjugate(PlusIOverlap))+"\n")
print("|<-i|psi>|^2 = "+str(MinusIOverlap*np.conjugate(MinusIOverlap))+"\n")

#aer_job = aer_sim.run(qc, shots=100)
#print(aer_job.result().get_counts())
#"""

qubitorder = [0,1,2,3,4,5,6]

KLZeroStates = []#[ZeroLogical]
KLOneStates = []#[OneLogical]

logicalStates = [ZeroLogical, OneLogical]
for j in range(2):
    for i in range(42):
        qc = qs.QuantumCircuit(7)

        qc.initialize(logicalStates[j], qc.qubits)

        if(i<21):
            if(i<7):
                qc.x(i)
            if(i>=7 and i < 14):
                qc.z(i%7)
            if(i>=14):
                qc.y(i%7)

        qc.unitary(mygate, qubitorder)

        if(i >= 21 and i < 42):
            if(i<28):
                qc.x(i%7)
            if(i>=28 and i < 35):
                qc.z(i%7)
            if(i>=35):
                qc.y(i%7)

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
    qc = qs.QuantumCircuit(7)
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