# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns


# def Calculate_lambda(A,U,V,m):
#     n = np.size(A,axis=0)
#     c = np.size(V,axis=0)
#     d = np.size(V,axis=1)

#     um = U**m

#     term1 = (um@V[:,:,0])/(um.sum(axis=1)[:,None] + 1e-8)
#     term1 = A[:,:,0] - term1

#     term1 = term1.sum(axis=0)   #check if axis=0
#     term2 = um.sum(axis=1)
#     term3 = (1/(term2 + 1e-8)).sum()

#     return 2*term1/term3



# def Update_B(A,V,U,m):
# 	um = U**m
#     term1 = (um@V[:,:,0])/(um.sum(axis=1)[:,None] + 1e-8)
#     term1 = A[:,:,0] - term1
#     term2 = Calculate_lambda(A,U,V,m)/(2*(um.sum(axis=1)[:,None] + 1e-8))
#     return term1-term2
    


# def Update_Y(A,B,g,beta=None,p=None):
#     n = np.size(A,axis=0)
#     d = np.size(A,axis=1)
#     Y = np.zeros((n,d,3))
#     Y[:,:,0] = A[:,:,0] - B
#     mu = Y[:,:,0]
#     if g == "yager":
#         Y[:,:,1] = (1 - mu**p)**(1/p)
#     elif g=="sugeno":
#         Y[:,:,1] = (1 - mu)/(1 + beta*mu)

#     Y[:,:,2] = 1-Y[:,:,0]-Y[:,:,1]

#     return Y
    

# def IFS_Distance_num_and_den(Y,V,i,j):
#     c = np.size(V,axis=0)
#     temp = 0.5 * np.sum((Y[j] - V[i])**2)

#     val = 0
#     for k in range(c):
#         temp2 = 0.5 * np.sum((Y[j] - V[k])**2)
#         val = val + temp2

#     return temp
    



# def Calculate_SNI(U,A,i,j,N,Y,V):
#     list1 = []
#     n = np.size(U,axis=0)
#     for k in range(n):
#         if(N[k][j] ==1):
#             list1.append(k)
#     val = 0
#     for t in list1:
#         val = val + U[i][t]

#     val2 = 0
#     for u in list1:
#         temp2 = IFS_Distance_num_and_den(Y,V,i,u)
#         val2 = val2 + (U[i][u]*temp2)/(val + 1e-8)
    
#     return val2




# def Update_U(A,V,B,a_,m,N,U,Y):
#     U_old = U
#     n = np.size(A,axis=0)
#     c = np.size(V,axis=0)
#     d = np.size(V,axis=1)
#     U = np.zeros((n,c))
#     for j in range(n):
#         for i in range(c):
#             temp = IFS_Distance_num_and_den(Y,V,i,j)
#             SNI = Calculate_SNI(U_old,A,i,j,N,Y,V)
#             den = 0
#             for t in range(c):
#                 temp = IFS_Distance_num_and_den(Y,V,t,j) + a_*Calculate_SNI(U_old,A,t,j,N,Y,V)
#                 den = den + temp**(-1/(m-1))
#             U[j][i] = ((temp + a_*SNI)**(-1/(m-1)))/den
            
            
# def Neighbourhood_Term_inside_membership(U,A,i,j,N,a_):
# 	n = np.size(U,axis=1)
#     list1 = []
#     for k_ in range(n):
#         if N[k_][j] == 1:
#             list1.append(k_)
#     den = 0
#     for k in list1:
#         den = den + U[i][k]

#     second_term = 0
#     for k in list1:
#         second_term = second_term + (U[i][k]*A[k,:,0])/den

#     return A[j,:,0] + a_*second_term
        
        
    
# def Neighbourhood_Term_inside_non_membership(U,A,i,j,N,a_):
# 	n = np.size(U,axis=1)
#     list1 = []
#     for k_ in range(n):
#         if N[k_][j] == 1:
#             list1.append(k_)
#     den = 0
#     for k in list1:
#         den = den + U[i][k]

#     second_term = 0
#     for k in list1:
#         second_term = second_term + (U[i][k]*A[k,:,1])/den

#     return A[j,:,1] + a_*second_term


# def Neighbourhood_Term_inside_hesitation(U,A,i,j,N,a_):
# 	n = np.size(U,axis=1)
#     list1 = []
#     for k_ in range(n):
#         if N[k_][j] == 1:
#             list1.append(k_)
#     den = 0
#     for k in list1:
#         den = den + U[i][k]

#     second_term = 0
#     for k in list1:
#         second_term = second_term + (U[i][k]*A[k,:,2])/den

#     return A[j,:,2] + a_*second_term


# def Update_V_membership(U,A,a_,m,N):
# 	n = np.size(A,axis=0)
#     d = np.size(A,axis=1)
#     c = np.size(U,axis=1)
#     V_membership = np.zeros((c,d))
#     um = U**m

#     for i in range(c):
#         num = 0
#         for j in range(n):
#             temp = (um[i][j])*(Neighbourhood_Term_inside_membership(U,A,i,j,N,a_))
#             num = num + temp
#         den = (1+a_)*um[i].sum()
#         V_membership[i] = num/den

#     return V_membership
    

    
    
    
# def Update_V_hesitation(U,A,a_,m,N):
# 	n = np.size(A,axis=0)
#     d = np.size(A,axis=1)
#     c = np.size(U,axis=1)
#     V_membership = np.zeros((c,d))
#     um = U**m

#     for i in range(c):
#         num = 0
#         for j in range(n):
#             temp = (um[i][j])*(Neighbourhood_Term_inside_hesitation(U,A,i,j,N,a_))
#             num = num + temp
#         den = (1+a_)*um[i].sum()
#         V_membership[i] = num/den

#     return V_membership



# def Update_V_non_membership(U,A,a_,m,N):
# 	n = np.size(A,axis=0)
#     d = np.size(A,axis=1)
#     c = np.size(U,axis=1)
#     V_membership = np.zeros((c,d))
#     um = U**m

#     for i in range(c):
#         num = 0
#         for j in range(n):
#             temp = (um[i][j])*(Neighbourhood_Term_inside_non_membership(U,A,i,j,N,a_))
#             num = num + temp
#         den = (1+a_)*um[i].sum()
#         V_membership[i] = num/den

#     return V_membership



# def Update_V(U,A,a_,m,N):
#     d = np.size(A,axis=1)
#     c = np.size(U,axis=1)
#     V = np.zeros((c,d,3))
#     V[:,:,0] = Update_V_membership(U,A,a_,m,N)
#     V[:,:,1] = Update_V_non_membership(U,A,a_,m,N)
#     V[:,:,2] = Update_V_hesitation(U,A,a_,m,N)
    
#     return V
    

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def Calculate_lambda(A,U,V,m):
    n = np.size(A,axis=0)
    c = np.size(V,axis=0)
    d = np.size(V,axis=1)

    um = U**m

    term1 = (um@V[:,:,0])/(um.sum(axis=1)[:,None] + 1e-8)
    term1 = A[:,:,0] - term1

    term1 = term1.sum(axis=0)   #check if axis=0
    term2 = um.sum(axis=1)
    term3 = (1/(term2 + 1e-8)).sum()

    return 2*term1.sum()/term3



def Update_B(A,V,U,m):
    um = U**m
    term1 = (um@V[:,:,0])/(um.sum(axis=1)[:,None] + 1e-8)
    term1 = A[:,:,0] - term1
    term2 = Calculate_lambda(A,U,V,m)/(2*(um.sum(axis=1)[:,None] + 1e-8))
    return term1-term2
    


def Update_Y(A,B,g,beta=None,p=None):
    n = np.size(A,axis=0)
    d = np.size(A,axis=1)
    Y = np.zeros((n,d,3))
    Y[:,:,0] = A[:,:,0] - B
    mu = Y[:,:,0]
    if g == "yager":
        Y[:,:,1] = (1 - mu**p)**(1/p)
    elif g=="sugeno":
        Y[:,:,1] = (1 - mu)/(1 + beta*mu)

    Y[:,:,2] = 1-Y[:,:,0]-Y[:,:,1]

    return Y
    


def IFS_Distance_num_and_den(Y,V,i,j):
    c = np.size(V,axis=0)
    temp = 0.5 * np.sum((Y[j] - V[i])**2)

    val = 0
    for k in range(c):
        temp2 = 0.5 * np.sum((Y[j] - V[k])**2)
        val = val + temp2

    return temp
    

def Calculate_SNI(U,A,i,j,N,Y,V):
    list1 = []
    n = np.size(U,axis=0)
    for k in range(n):
        if(N[k][j] ==1):
            list1.append(k)
    val = 0
    for t in list1:
        val = val + U[i][t]

    val2 = 0
    for u in list1:
        temp2 = IFS_Distance_num_and_den(Y,V,i,u)
        val2 = val2 + (U[i][u]*temp2)/(val + 1e-8)
    
    return val2



def Update_U(A,V,B,a_,m,N,U,Y):
    U_old = U
    n = np.size(A,axis=0)
    c = np.size(V,axis=0)
    d = np.size(V,axis=1)
    U = np.zeros((n,c))
    for j in range(n):
        for i in range(c):
            temp = IFS_Distance_num_and_den(Y,V,i,j)
            SNI = Calculate_SNI(U_old,A,i,j,N,Y,V)
            den = 0
            for t in range(c):
                temp = IFS_Distance_num_and_den(Y,V,t,j) + a_*Calculate_SNI(U_old,A,t,j,N,Y,V)
                den = den + temp**(-1/(m-1))
            U[j][i] = ((temp + a_*SNI)**(-1/(m-1)))/den

    return U


            
            

def Neighbourhood_Term_inside_membership(U,A,i,j,N,a_):
	n = np.size(U,axis=1)
    list1 = []
    for k_ in range(n):
        if N[k_][j] == 1:
            list1.append(k_)
    den = 0
    for k in list1:
        den = den + U[i][k]

    second_term = 0
    for k in list1:
        second_term = second_term + (U[i][k]*A[k,:,0])/den

    return A[j,:,0] + a_*second_term
        
        
    
def Neighbourhood_Term_inside_non_membership(U,A,i,j,N,a_):
	n = np.size(U,axis=1)
    list1 = []
    for k_ in range(n):
        if N[k_][j] == 1:
            list1.append(k_)
    den = 0
    for k in list1:
        den = den + U[i][k]

    second_term = 0
    for k in list1:
        second_term = second_term + (U[i][k]*A[k,:,1])/den

    return A[j,:,1] + a_*second_term


def Neighbourhood_Term_inside_hesitation(U,A,i,j,N,a_):
	n = np.size(U,axis=1)
    list1 = []
    for k_ in range(n):
        if N[k_][j] == 1:
            list1.append(k_)
    den = 0
    for k in list1:
        den = den + U[i][k]

    second_term = 0
    for k in list1:
        second_term = second_term + (U[i][k]*A[k,:,2])/den

    return A[j,:,2] + a_*second_term



def Update_V_membership(U,A,a_,m,N):
    d = np.size(A,axis=1)
    c = np.size(U,axis=1)
    V_membership = np.zeros((c,d))
    um = U**m

    for i in range(c):
        num = 0
        for j in range(n):
            temp = (um[j][i])*(Neighbourhood_Term_inside_membership(U,A,i,j,N,a_))
            num = num + temp
        den = (1+a_)*um[i].sum()
        V_membership[i] = num/den

    return V_membership
    

def Update_V_hesitation(U,A,a_,m,N):
    d = np.size(A,axis=1)
    c = np.size(U,axis=1)
    V_membership = np.zeros((c,d))
    um = U**m

    for i in range(c):
        num = 0
        for j in range(n):
            temp = (um[j][i])*(Neighbourhood_Term_inside_hesitation(U,A,i,j,N,a_))
            num = num + temp
        den = (1+a_)*um[i].sum()
        V_membership[i] = num/den

    return V_membership
    
    

def Update_V_non_membership(U,A,a_,m,N):
    d = np.size(A,axis=1)
    c = np.size(U,axis=1)
    V_membership = np.zeros((c,d))
    um = U**m

    for i in range(c):
        num = 0
        for j in range(n):
            temp = (um[j][i])*(Neighbourhood_Term_inside_non_membership(U,A,i,j,N,a_))
            num = num + temp
        den = (1+a_)*um[i].sum()
        V_membership[i] = num/den

    return V_membership
    



def Update_V(U,A,a_,m,N):
    d = np.size(A,axis=1)
    c = np.size(U,axis=1)
    V = np.zeros((c,d,3))
    V[:,:,0] = Update_V_membership(U,A,a_,m,N)
    V[:,:,1] = Update_V_non_membership(U,A,a_,m,N)
    V[:,:,2] = Update_V_hesitation(U,A,a_,m,N)

    return V
    
    


def Calculate_D_ifs(muA,nuA,piA,muV,nuV,piV,i,j):

    term1 = np.linalg.norm(muA[j] - muV[i])**2
    term2 = np.linalg.norm(nuA[j] - nuV[i])**2
    term3 = np.linalg.norm(piA[j] - piV[i])**2

    return term1 + term2 + term3


U_old = U.copy()
V_old = V.copy()
B_old = B.copy()
Y_old = Y.copy()



itn = 0

while(itn < 100):
    Y_new = Update_Y(A,B_old,g,beta=None,p=None)

    U_new = Update_U(A,V_old,B_old,a_,m,N,U_old,Y_new)

    V_new = Update_V(U_new,A,a_,m,N)

    B_new = Update_B(A,V_new,U_new,m)

    if Diff(U_old,U_new,V_old,V_new,B_old,B_new) < e:
        break
   
    U_old = U_new.copy()
    B_old = B_new.copy()
    V_old = V_new.copy()
    Y_old = Y_new.copy()

    itn = itn + 1