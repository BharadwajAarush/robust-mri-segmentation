import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

c = int(input("Enter the number of clusters"))
n = int(input("Enter the number of points"))
d = int(input("Enter the dimension of points"))
m = float(input("Enter the value of the fuzzifier (recommended >1)"))
e = float(input("Enter the Value of Epsilon"))
a_ = float(input("Enter the value of the spatial constant"))


A = np.random.rand(n,d)
V = np.random.rand(c,d)
B = np.random.rand(n,d)
U = np.random.rand(n,c)
N = np.random.randint(2,size(n,n))
N = N + N.T
np.fill_diagonal(N, 0)
for i in range(n):
    sum_i = U[i].sum()
    U[i] = U[i]/sum_i



def Update_V(A,B,U,m,a_,N):
    um = U**m
    term1 = um.T@(A-B)
    c = np.size(U,axis=1)
    d = np.size(A,axis=1)
    V = np.zeros((c,d))

    for i in range(c):
        V[i] = (term1 + term2_in_V(A,B,U,m,a_,N,i))/(1+a_)*um.sum(axis=0)[:,None]

    return V
    
    
    

def term2_in_V(A,B,U,m,a_,N,i):
    xr = np.empty((1,d))
    br = np.empty((1,d))
    um = U**m
    for j in range(n):
        if N[i][j]==1:
            np.append(xr,A[j],axis=0)

    for j in range(n):
        if N[i][j] == 1:
            np.append(br,B[j],axis=0)


    nr = np.size(br,axis=0)

    term2 = (a_/nr)*(um.T@(xr-br))
    
            


rough = np.array([[1,2,3],[4,5,6]])
print(np.size(rough,axis=0))# no of rows
print(np.size(rough,axis=1))#no of columns/features


def Calculate_Sum_for_U(A,B,V,m,a_,N,i,j):
    term1 = np.linalg.norm(A[j]-B[j]-V[i])
    term1 = term1**2
    list1 = []
    for r in range(n):
        if N[r][j] ==1:
            list.append(r)

    val = 0;
    nr = 0

    for t in list1:
        nr = nr+1
        temp = np.linalg.norm(A[t]-B[t] -V[i])
        val = val + temp**2

    term2 = (a_/nr)*val

    return term1+term2
        
            
    
def Calculate_Denominator_for_U(A,B,V,m,a_,N,j):
    c = np.size(V,axis=0)
    val = 0
    for i in range(c):
        temp = Calculate_Sum_for_U(A,B,V,m,a_,N,i,j)
        val = val + temp**(-1/m-1)
    
        
    


def Update_U(A,B,V,m,a_,N):
    n = np.size(A,axis=0)
    c = np.size(V,axis=0)

    U = np.eros((n,c))
    for i in range(n):
        for j in range(c):
            num = Calculate_Sum_for_U(A,B,V,m,a_,N,i,j)
            den = Calculate_Denominator_for_U(A,B,V,m,a_,N,j):
            U[i][j] = (num**(-1/m-1))/den
    
    
def Update_B(A,V,U,m):
    um = U**m
    B = A-((um@V)/(um.sum(axis=1)[:,None]))
    return B

## done according to the paper , not according to reality 

def Diff(U_old,U_new,V_old,V_new,B_old,B_new):
    val1 = np.linalg.norm(U_old-U_new,2)
    val2 = np.linalg.norm(V_new-V_old,2)
    val3 = np.linalg.norm(B_new-B_old,2)
    return val1+val2+val3



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
