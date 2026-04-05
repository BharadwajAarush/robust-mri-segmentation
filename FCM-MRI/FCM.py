# THOUGHT BEHIND THE IMPLEMENTATION-
# 1- WE WILL BUILD A MEMBERSHIP MATRIX CALLED U
# 2- WE WILL BUILD A PARTITION MATRIX 
# 3- WE WILL BUILD A MATRIX OF ALL THE CENTERS CALLED V


# ALGORITHM-
# 1- TAKE THE VALUES OF M(FUZZIFIER > 1) , C(NUMBER OF CLUSTERS ) , E(EPSILON , CORRESONDING TO THE MINIMUM ERROR ALLOWED FOR CONVERGENCE) 
# 2- START WITH A MATRIX 'A' OF DIMENSION nXm SUCH THAT n CORRESPONDS TO THE NUMBER OF POINTS IN THE ARRAY AND m CORRESPONDS TO THE DIFFERENT FEATURES OF IT ( USUALLY IT IS DIFF VALUES OF RGB FOR INTENSITY , OR DIFFERENT VALUES OF x, y,z COORDINATES FOR SPATIAL ALLIGNMENT
# 3- START WITH A MATRIX 'V' of dimension cXm , IT CONTAINS ALL THE CENTER POINTS OF THE C CLUSTERS , START WITH A RANDOM SET OF POINTS
# 4- START WITH A MATRIX 'U' of dimension n*c , IT CONTAINS THE MEMBERSHIP VALUE OF EACH POINT TO THE CORRESPONSDING CLUSTER CENTER , START WITH A RANDOM U MATRIX SUCH THAT , SUMMATION OF ALL THE VALUES IN A ROW IS 1
# 5- FOR A FIXED U , UPDATE V WITH A CERTAIN FORMULA
# 6- FOR A FIXED V , UPDATE U WITH A CERTAIN FORMULA 
# 7- REPEAT THE STEPS TILL WE CONVERGE , OR REPEAT TILL THERE IS NO CHANGE IN THE MATRIX U , OR NORM OF CHANGE BETWEEN THE OLD U AND NEW U IS MINIMAL 
# 8- WITH THE SPOT ON SET OF U AND V , MAKE ANOTHER MATRIX CALLED 'PARTITION' , AND PARTITION ALL THE POINTS , SOFTLY BASED ON THE NEW MATRICES
# FUNCTIONS TO DEFINE-
# 1- UPDATE U , the inputs of the given function are updated value of V and m 
# 2- UPDATE V , the imputs of the given function are updated values of U and m and c
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
def update_V(U,A,m):
    U = np.array(U)
    um = U**m
    V = (um.T@A)/um.sum(axis=0)[:,None]
    return V
    
# def update_U(A,V,m):
#     n = np.size(A,axis=0)
#     c = np.size(V,axis=0)
#     for i in range(n):
#         val = 0
#         for j in range(c):
#             val = val + np.pow((np.linalg.norm(A[i] - V[j])),-2/(m-1))
#         for k in range(c):
#             U[i][k] =  np.pow((np.linalg.norm(A[i] - V[j])),-2/(m-1)/)val
#     return U

def update_U(A,V,m):
    n = np.size(A,axis=0)
    c = np.size(V,axis=0)
    U = np.zeros((n,c))
    for i in range(n):
        val = 0
        for j in range(c):
            temp = np.linalg.norm(A[i]-V[j])**2
            val = val+ temp**(-2/(m-1))
            # val = val + np.pow((np.linalg.norm(A[i] - V[j]),2),-2/(m-1))
        for k in range(c):
            temp2 = np.linalg.norm(A[i]- V[k])**2
            U[i,k] = temp2**(-2/(m-1))/val
            # U[i,k] =  np.pow((np.linalg.norm(A[i] - V[k]),2),-2/(m-1))/val
    return U

c = int(input("Enter the number of clusters"))
n = int(input("Enter the number of points"))
d = int(input("Enter the dimension of points"))
m = float(input("Enter the value of the fuzzifier (recommended >1)"))
e = float(input("Enter the Value of Epsilon"))


# def Diff(U_old,U_new,V_old,V_new):
#     val1 = np.linalg.norm(U_old-U_new)
#     val2 = np.linalg.norm(V_new-V_old)
#     return val1+val2

def Diff(U_old,U_new,V_old,V_new):
    val1 = np.linalg.norm(U_old-U_new,2)
    val2 = np.linalg.norm(V_new-V_old,2)
    return val1+val2


A = np.random.rand(n,d)
V = np.random.rand(c,d)
U = np.random.rand(n,c)
for i in range(n):
    sum_i = U[i].sum()
    U[i] = U[i]/sum_i


# U_old = U
# V_old = V
# U_new = update_U(A,V_old,m)
# V_new = update_V(U_new,A,m)
# while(Diff(U_old,U_new,V_old,V_new)>e):
#     U_old = U_new
#     V_old = V_new
#     U_new = update_U(A,V_old,m)
#     V_new = update_V(U_new,A,m)


U_old = U.copy()
V_old = V.copy()
U_new = update_U(A,V_old,m)
V_new = update_V(U_new,A,m)
itn = 0
while(Diff(U_old,U_new,V_old,V_new)>e AND itn>1000):
    U_old = U_new.copy()
    V_old = V_new.copy()
    U_new = update_U(A,V_old,m)
    V_new = update_V(U_new,A,m)
    itn = itn +1