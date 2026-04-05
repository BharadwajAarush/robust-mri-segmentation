import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def Update_v(A,B,U,m):
    um = U**m
    V = um.T@(A-B)/um.sum(axis=0)[:,None]
    return V


def Update_B(A,V,U,m):
    um = U**m
    B = A-((um@V)/(um.sum(axis=1)[;,None]))
    return B


def Update_U(A,B,V,m):
    n = np.size(A,axis=0)
    c = np.size(V,axis=0)
    d = np.size(A,axis=1)
    U = np.zeros((n,c))
    for i in range(n):
        val=0
        for j in range(c):
            temp = np.linalg(A[i]-B[i]-V[j],2)
            val = val + temp**((-1)/(m-1))
        for k in range(c):
            temp2 = np.linalg(A[i]-B[i]-V[k],2)
            U[i][k] = (temp2**(-1/(m-1)))/val

    return U



c = int(input("Enter the number of clusters"))
n = int(input("Enter the number of points"))
d = int(input("Enter the dimension of points"))
m = float(input("Enter the value of the fuzzifier (recommended >1)"))
e = float(input("Enter the Value of Epsilon"))



def Diff(U_old,U_new,V_old,V_new,B_old,B_new):
    val1 = np.linalg.norm(U_old-U_new,2)
    val2 = np.linalg.norm(V_new-V_old,2)
    val3 = np.linalg.norm(B_new-B_old,2)
    return val1+val2+val3



A = np.random.rand(n,d)
V = np.random.rand(c,d)
B = np.random.rand(n,d)
U = np.random.rand(n,c)




for i in range(n):
    sum_i = U[i].sum()
    U[i] = U[i]/sum_i



U_old = U.copy()
V_old = V.copy()
B_old = B.copy()
V_new = update_V(A,B_old,U_old,m)
U_new = update_U(A,B_old,V_new,m)
B_new = update_B(A,V_new,U_new,m)
itn = 0


while(Diff(U_old,U_new,V_old,V_new)>e and itn>1000):
    U_old = U_new.copy()
    V_old = V_new.copy()
    B_old = B_new.copy()
    V_new = update_V(A,B_old,U_old,m)
    U_new = update_U(A,B_old,V_new,m)
    B_new = update_B(A,V_new,U_new,m)
    
    
    itn = itn + 1
