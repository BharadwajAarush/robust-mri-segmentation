# A LOT LIKE FCM WITH ADDITIONAL 'B' WHICH CAN BE CATERED DOWN WITH EASE
# Algorithm(Alternate Minimization Stratergy)-
# 1-UPDATE V , FIXING U AND B
# 2-UPDATE U , FIXING V AND B
# 3-UPDATE B , FIXING V AND U

# ALOT LIKE BCFCM WITHOUT SPATIAL BUT WE HAVE TO FORM THIS NEIGHBOURHOOD MATRIX AS WELL
# THE PROPERTIES OF NEIGHBOURHOOD MATRIX MUST BE -
# 1- IS A POINT IS NEIGHBOUR OF OTHER THEN THE OTHER POINT MUST BE IT'S NEIGHBOUR AS WELL
# 2- TO RANDOMLY GENERATE IT , MAKE A SYMMETRIC MATRIX OF ALL NON ZERO ENTRIES 1;
# WHERE 1 SYMBOLISES NEIGHBOUR AND 0 SYMBOLIZES NOT NEIGHBOUR

import numpy as np
from scipy.signal import convolve2d


def Build_Neighborhood_Kernel(window=3):
    return np.ones((window, window)) / (window * window)


def Update_V(U, A, B, x_bar, a_, m):
    um = U ** m
    # term1 + term2_in_V
    x_eff = (A - B) + a_ * (x_bar - B)
    V = (um @ x_eff) / ((1.0 + a_) * um.sum(axis=1)[:, None] + 1e-8)
    return V


def Update_U(A, B, x_bar, V, a_, m):
    # Calculate_Sum_for_U
    diff1 = (A - B).T - V   # (c, n)
    diff2 = (x_bar - B).T - V

    cost = diff1 ** 2 + a_ * (diff2 ** 2)
    cost = np.clip(cost, 1e-10, None)

    # Calculate_Denominator_for_U
    expo = -1.0 / (m - 1.0)
    cost_pow = cost ** expo
    den = cost_pow.sum(axis=0, keepdims=True)
    U = cost_pow / (den + 1e-8)
    return U


def Update_B(A, V, U, x_bar, a_, m, h, w, bias_kernel):
    um = U ** m
    v_eff = (um.T @ V) / (um.sum(axis=0)[:, None] + 1e-8)
    raw_B = (A + a_ * x_bar) / (1.0 + a_) - v_eff

    ## done according to the paper , not according to reality 
    # smooth the bias field with a 5x5 window because RF bias field is physically smooth
    B_2d = raw_B.reshape(h, w)
    B_smooth = convolve2d(B_2d, bias_kernel, mode='same', boundary='symm')
    return B_smooth.reshape(-1, 1)


def Diff(U_old, U_new, V_old, V_new, B_old, B_new):
    val1 = np.linalg.norm(U_new - U_old)
    val2 = np.linalg.norm(V_new - V_old)
    val3 = np.linalg.norm(B_new - B_old)
    return val1 + val2 + val3


def Run(img, c=3, m=2.0, a_=1.0, e=1e-4, max_iter=50, seed=42):
    np.random.seed(seed)
    img_norm = img.astype(float)
    img_norm = (img_norm - img_norm.min()) / (img_norm.max() - img_norm.min() + 1e-8)

    h, w = img.shape
    n = h * w
    A = img_norm.reshape(n, 1)

    k3 = Build_Neighborhood_Kernel(3)
    k5 = Build_Neighborhood_Kernel(5)
    x_bar = convolve2d(img_norm, k3, mode='same', boundary='symm').reshape(n, 1)

    quantiles = np.linspace(0.1, 0.9, c)
    V = np.quantile(A, quantiles).reshape(c, 1)
    B = np.zeros((n, 1))

    U = Update_U(A, B, x_bar, V, a_, m)

    U_old = U.copy()
    V_old = V.copy()
    B_old = B.copy()
    itn = 0

    while itn < max_iter:
        V_new = Update_V(U_old, A, B_old, x_bar, a_, m)
        U_new = Update_U(A, B_old, x_bar, V_new, a_, m)
        B_new = Update_B(A, V_new, U_new, x_bar, a_, m, h, w, k5)

        if Diff(U_old, U_new, V_old, V_new, B_old, B_new) < e:
            U_old = U_new
            V_old = V_new
            B_old = B_new
            break

        U_old = U_new.copy()
        V_old = V_new.copy()
        B_old = B_new.copy()
        itn = itn + 1

    labels = np.argmax(U_old, axis=0)
    seg = labels.reshape(h, w)

    return seg, U_old, V_old, B_old
