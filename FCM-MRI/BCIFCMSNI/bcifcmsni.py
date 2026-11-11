# BCIFCMSNI (Bias-Corrected Intuitionistic Fuzzy C-Means with Spatial Neighborhood Information)
# Implementation based on BC-IFS-SNI-VANILLA.py and bcifswsni_bn.py

import numpy as np
from scipy.signal import convolve2d


def Image_to_IFS(img, g="yager", p=1.0, beta=1.5):
    img_norm = img.astype(float)
    img_norm = (img_norm - img_norm.min()) / (img_norm.max() - img_norm.min() + 1e-8)

    n = img_norm.shape[0] * img_norm.shape[1]
    d = 1

    mu = img_norm.reshape(n, d)
    mu = np.clip(mu, 0.0, 1.0)

    if g == "yager":
        nu = (1.0 - mu ** p) ** (1.0 / p)
    elif g == "sugeno":
        nu = (1.0 - mu) / (1.0 + beta * mu + 1e-8)
    else:
        nu = 1.0 - mu

    nu = np.clip(nu, 0.0, 1.0)
    pi = np.clip(1.0 - mu - nu, 0.0, 1.0)

    A = np.zeros((n, d, 3))
    A[:, :, 0] = mu
    A[:, :, 1] = nu
    A[:, :, 2] = pi
    return A


def Update_Y(A, B, g="yager", p=1.0, beta=1.5):
    n = np.size(A, axis=0)
    d = np.size(A, axis=1)
    Y = np.zeros((n, d, 3))
    Y[:, :, 0] = np.clip(A[:, :, 0] - B, 0.0, 1.0)
    mu_Y = Y[:, :, 0]

    if g == "yager":
        Y[:, :, 1] = (1.0 - mu_Y ** p) ** (1.0 / p)
    elif g == "sugeno":
        Y[:, :, 1] = (1.0 - mu_Y) / (1.0 + beta * mu_Y + 1e-8)
    else:
        Y[:, :, 1] = 1.0 - mu_Y

    Y[:, :, 1] = np.clip(Y[:, :, 1], 0.0, 1.0)
    Y[:, :, 2] = np.clip(1.0 - Y[:, :, 0] - Y[:, :, 1], 0.0, 1.0)
    return Y


def Calculate_D_ifs(Y, V):
    # IFS_Distance_num_and_den / Calculate_D_ifs
    c = V.shape[0]
    n = Y.shape[0]
    D = np.zeros((c, n))
    for i in range(c):
        diff = Y - V[i]
        D[i] = 0.5 * np.sum(diff ** 2, axis=(1, 2))
    return D


def Calculate_SNI(U, D, h, w, kernel):
    # Calculate_SNI: spatial neighborhood information weighted by U
    c, n = U.shape
    SNI = np.zeros((c, n))
    for i in range(c):
        U2d = U[i].reshape(h, w)
        den2d = convolve2d(U2d, kernel, mode='same', boundary='fill', fillvalue=0)
        num2d = convolve2d((U[i] * D[i]).reshape(h, w), kernel, mode='same', boundary='fill', fillvalue=0)
        SNI[i] = (num2d / (den2d + 1e-8)).ravel()
    return SNI


def Update_U(A, V, B, a_, m, U, Y, h, w, kernel):
    D = Calculate_D_ifs(Y, V)
    SNI = Calculate_SNI(U, D, h, w, kernel)

    cost = D + a_ * SNI
    cost = np.clip(cost, 1e-10, None)

    expo = -1.0 / (m - 1.0)
    cost_pow = cost ** expo
    den = cost_pow.sum(axis=0, keepdims=True)
    U_new = cost_pow / (den + 1e-8)
    return U_new


def Update_V(U, A, a_, m, h, w, kernel):
    # Update_V: membership, non-membership, hesitation
    c = U.shape[0]
    n = h * w
    um = U ** m
    V = np.zeros((c, 1, 3))

    for comp in range(3):
        A2d = A[:, 0, comp].reshape(h, w)
        for i in range(c):
            U2d = U[i].reshape(h, w)
            den2d = convolve2d(U2d, kernel, mode='same', boundary='fill', fillvalue=0)
            num2d = convolve2d(U2d * A2d, kernel, mode='same', boundary='fill', fillvalue=0)
            second2d = num2d / (den2d + 1e-8)
            eff2d = A2d + a_ * second2d
            num = (um[i].reshape(h, w) * eff2d).sum()
            den_tot = (1.0 + a_) * um[i].sum()
            V[i, 0, comp] = num / (den_tot + 1e-8)
    return V


def Update_B(A, V, U, h, w, bias_kernel):
    # term1 = term1.sum(axis=0)   #check if axis=0
    um = U ** 2.0
    v_eff = (um.T @ V[:, :, 0]) / (um.sum(axis=0)[:, None] + 1e-8)
    raw_B = A[:, :, 0] - v_eff

    # smooth bias field over 5x5 window
    B_2d = raw_B.reshape(h, w)
    B_smooth = convolve2d(B_2d, bias_kernel, mode='same', boundary='symm')
    return B_smooth.reshape(-1, 1)


def Diff(U_old, U_new, V_old, V_new, B_old, B_new):
    val1 = np.linalg.norm(U_new - U_old)
    val2 = np.linalg.norm(V_new - V_old)
    val3 = np.linalg.norm(B_new - B_old)
    return val1 + val2 + val3


def Run(img, c=4, m=2.0, a_=3.0, g="yager", p=1.0, beta=1.5, e=1e-4, max_iter=50, seed=42):
    np.random.seed(seed)
    A = Image_to_IFS(img, g=g, p=p, beta=beta)

    h, w = img.shape
    n = h * w

    quantiles = np.linspace(0.1, 0.9, c)
    mu_centers = np.quantile(A[:, 0, 0], quantiles)

    V = np.zeros((c, 1, 3))
    for i in range(c):
        mu_val = mu_centers[i]
        if g == "yager":
            nu_val = (1.0 - mu_val ** p) ** (1.0 / p)
        elif g == "sugeno":
            nu_val = (1.0 - mu_val) / (1.0 + beta * mu_val + 1e-8)
        else:
            nu_val = 1.0 - mu_val
        pi_val = max(0.0, 1.0 - mu_val - nu_val)
        V[i, 0, 0] = mu_val
        V[i, 0, 1] = nu_val
        V[i, 0, 2] = pi_val

    B = np.zeros((n, 1))  # small values also fine
    Y = Update_Y(A, B, g=g, p=p, beta=beta)

    D = Calculate_D_ifs(Y, V)
    D = np.clip(D, 1e-10, None)
    cost_init = D ** (-1.0 / (m - 1.0))
    U = cost_init / (cost_init.sum(axis=0, keepdims=True) + 1e-8)

    k3 = np.ones((3, 3))
    np.fill_diagonal(k3, 1.0)
    k5 = np.ones((5, 5)) / 25.0

    U_old = U.copy()
    V_old = V.copy()
    B_old = B.copy()
    itn = 0

    while itn < max_iter:
        Y_new = Update_Y(A, B_old, g=g, p=p, beta=beta)
        U_new = Update_U(A, V_old, B_old, a_, m, U_old, Y_new, h, w, k3)
        V_new = Update_V(U_new, A, a_, m, h, w, k3)
        B_new = Update_B(A, V_new, U_new, h, w, k5)

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
