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


def Update_V(U, A, m):
    um = U ** m
    V = (um @ A) / (um.sum(axis=1)[:, None] + 1e-8)
    return V


def Update_U(A, V, m):
    diff = A[:, None, :] - V[None, :, :]
    dist = np.sum(diff ** 2, axis=2).T
    dist = np.clip(dist, 1e-10, None)

    expo = -1.0 / (m - 1.0)
    cost = dist ** expo
    den = cost.sum(axis=0, keepdims=True)
    U = cost / (den + 1e-8)
    return U


def Diff(U_old, U_new, V_old, V_new):
    val1 = np.linalg.norm(U_new - U_old)
    val2 = np.linalg.norm(V_new - V_old)
    return val1 + val2


def Run(img, c=3, m=2.0, e=1e-4, max_iter=50, seed=42):
    np.random.seed(seed)
    img_norm = img.astype(float)
    img_norm = (img_norm - img_norm.min()) / (img_norm.max() - img_norm.min() + 1e-8)

    h, w = img.shape
    n = h * w
    A = img_norm.reshape(n, 1)

    quantiles = np.linspace(0.1, 0.9, c)
    V = np.quantile(A, quantiles).reshape(c, 1)

    U = Update_U(A, V, m)

    U_old = U.copy()
    V_old = V.copy()
    itn = 0

    while itn < max_iter:
        V_new = Update_V(U_old, A, m)
        U_new = Update_U(A, V_new, m)

        if Diff(U_old, U_new, V_old, V_new) < e:
            U_old = U_new
            V_old = V_new
            break

        U_old = U_new.copy()
        V_old = V_new.copy()
        itn = itn + 1

    labels = np.argmax(U_old, axis=0)
    seg = labels.reshape(h, w)

    return seg, U_old, V_old
