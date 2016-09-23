## @package "finite-volume-method"
#  This is a package for the computation of 1D finite volume method.
#
#  Two techniques are implemented: simple upwind scheme and weighted
#  essentially non-oscillatory (WENO) scheme.

import numpy as np
import sys
import re
import math

cellWidth = 2

## Simple upwind scheme.
#
# More details.
# Some in-text formula: \f$ \sigma=\int_V \partial x/\partial x \f$
# Some new-line formula
# \f[
#   |I_2|=
#     \left| \int_{0}^T \psi(t) 
#     \left\{ 
#     u(a,t)-
#     \int_{\gamma(t)}^a 
#     \frac{d\theta}{k(\theta,t)}
#     \int_{a}^\theta c(\xi)u_t(\xi,t)\,d\xi
#     \right\} dt
#     \right|
# \f]
# \f{eqnarray*}{
#     g &=& \frac{Gm_2}{r^2} \\ 
#     &=& \frac{(6.673 \times 10^{-11}\,\mbox{m}^3\,\mbox{kg}^{-1}\,
#                \mbox{s}^{-2})(5.9736 \times 10^{24}\,\mbox{kg})}{(6371.01\,\mbox{km})^2} \\ 
#     &=& 9.82066032\,\mbox{m/s}^2
# \f}
def advect(u,v,dt,dx):

    numCells = len(u)
    u_right = np.linspace(0,0,numCells)
    u_left = np.copy(u_right)

    u_out = np.copy(u)

    for i in range(len(v)):

        ip = i+1 if i < len(v)-1 else i
        im = i-1 if i > 0 else 0

        v_right = (v[i]+v[ip])/2.0
        v_left = (v[i]+v[im])/2.0

        right = 0.0 if v_right < 0.0 else 1.0
        left = 0.0 if v_left > 0.0 else 1.0

        u_right[i] =  v_right * u[i] * right;
        u_left[i] = v_left * u[i] * left;

    residual = np.linspace(0,0,numCells+1)

    for i in range(len(u)):
        flux = u_right[i] * dt / dx
        residual[i+1] += flux
        residual[i] -= flux

    for i in range(len(u)):
        u_out[i] += residual[i]
        
    residual = np.linspace(0,0,numCells+1)

    for i in range(len(u)):
        flux = - u_left[i] * dt / dx
        residual[i] += flux
        residual[i+1] -= flux

    for i in range(len(u)):
        u_out[i] += residual[i+1]
        
    return u_out

def generateStencilMat(polyOrder):
    ### it is assumed that one cell has size of 2 units ###
    
    # L is the stencil matrix that can be procomputed
    L = np.ndarray(shape=(polyOrder+1,polyOrder+1,polyOrder+1))

    # 1,1/2,1/3 - come from the polynom integral (\int ax^2 = a/3 x^3, and so on)
    coeffs = np.array([1.0/x for x in range(1,polyOrder+2)])

    # intervals of integration
    intervals = np.ndarray(shape=(polyOrder+1,2))

    for k in range(polyOrder+1):
        
        cell = int(-polyOrder/2+k) #!! what about odd degrees?
        
        # if polynomial is of even order, then only one stencil is needed
        if polyOrder % 2 == 0:
            for i in range(-int(polyOrder/2),int(polyOrder/2)+1):
                intervals[i+int(polyOrder/2)][0] = ((i+cell)*2-1)
                intervals[i+int(polyOrder/2)][1] = ((i+cell)*2+1)
            
        for i in range(polyOrder+1):
            for j in range(polyOrder+1):
                L[k][i][j] = coeffs[j]*(intervals[i][1]**(j+1) - intervals[i][0]**(j+1))

    return L

## Oscillation matrix.
#
# \f[
#   \Sigma_{mn}=\sum_{r=1}^p\int\frac{\partial^r x^m}{\partial x^r}\frac{\partial^r x^n}{\partial x^r} \mathrm{d} x
# \f]
def generateOscillMat(polyOrder, dx):

    S = np.zeros(shape=(polyOrder,polyOrder))
    oscillMatSize = polyOrder

    for m in range(1,polyOrder+1):
        for n in range(1,polyOrder+1):
            s_mn = 0.0
            for r in range(1,polyOrder+1):

                scale = (1.0/dx)**r # 2/dx or 1/dx ??
                
                if r <= m and r <= n:
                    x_exp = m-r + n-r # exponent after multiplying derivatives
                    x_coeff = 1.0 # coefficient after multiplying derivatives
                    for c in range(r):
                        x_coeff *= (m-c)*(n-c)
                    int_x_exp = x_exp+1
                    int_x_coeff = x_coeff/int_x_exp

                    integral = int_x_coeff*1**int_x_exp - int_x_coeff*(-1)**int_x_exp
                    s_mn += integral * scale
                    
            S[m-1][n-1] = s_mn

    return S

def calcOmega(polyCoeffs, S, polyOrder):
    # polyCoeffs: polynomial coefficients; shape(polyCoeffs) = (polyOrder+1,polyOrder+1)
    # S: oscillation matrix; shape(S) = (polyOrder,polyOrder)
    
    omega = np.linspace(0.0,0.0,polyOrder+1)
    sumOmega = 0.0

    for j in range(polyOrder+1):

        states = S.dot(polyCoeffs[j][1:])
        sigma = states.dot(polyCoeffs[j][1:])

        la = 1000 if j == int(polyOrder/2) else 1 #!! what about odd degrees?
        omega[j] = la/((1.e-5 + sigma)**4)
        sumOmega += omega[j]

    if sumOmega > 0:
        omega = omega/sumOmega

    return omega

# Generates polynomial coefficients.

# More details.
# \f[
# \f{pmatrix}
# \hat{w}_0^{k} \\
# \hat{w}_1^{k} \\
# \hat{w}_2^{k}
# \f}
# = \mathbf{L}_{k}^{-1}
# \f{pmatrix}
# \overline{\phi}_{i-1}\\
# \overline{\phi}_{i}\\
# \overline{\phi}_{i+1}
# \f}
# \f]
def calcPolyCoeffs(polyOrder, f, i, Linv):

    # we take into account polyOrder+1 polynomials, each of them with
    # polyOrder+1 coefficients
    
    numCells = len(f)    
    polyCoeffs = np.ndarray(shape=(polyOrder+1,polyOrder+1))
    for j in range(polyOrder+1):

        # ids is used again to avoid if-else blocks regarding values on
        # the boundary; it clamps the indices to [0:numCells-1]
        ids = list(range(i-polyOrder+j,i+j+1))
        for k in range(len(ids)):
            ids[k] = max(0, min(ids[k], numCells-1))

        scale = cellWidth
        phi = np.array([scale*f[x] for x in ids])
        polyCoeffs[j] = Linv[j].dot(phi)

    return polyCoeffs

def advectWENO(f,u,dt,dx,polyOrder):

    numCells = len(f)
    
    f_right = np.linspace(0,0,numCells)
    f_left = np.copy(f_right)

    f_out = np.copy(f)

    S = generateOscillMat(polyOrder, dx)
    Linv = np.linalg.inv(generateStencilMat(polyOrder))

    for i in range(len(u)):

        # clamp velocity on the boundary
        ip = i+1 if i < len(u)-1 else i
        im = i-1 if i > 0 else 0

        # interpolate velocity at face between cells i and i+1
        u_right = (u[i]+u[ip])/2.0
        # interpolate velocity at face between cells i-1 and i
        u_left = (u[i]+u[im])/2.0

        #
        # only take into account fluxes out of the cell
        # this is done becuase we know the WENO reconstruction only at the
        # current cell and we use it in the upwind scheme (i.e., the current
        # cell is the "upwind" one only if velocity goes out of it)
        # to avoid mutliple if-else blocks, each flux is multiplied by the
        # corresponding left/rigth flag
        # 
        right = 0.0 if u_right < 0.0 else 1.0
        left = 0.0 if u_left > 0.0 else 1.0

        polyCoeffs = calcPolyCoeffs(polyOrder, f, i, Linv)
        omegas = calcOmega(polyCoeffs, S, polyOrder)
        wenoCoeffs = np.transpose(polyCoeffs).dot(omegas)

        for k in range(0,polyOrder+1):
            f_right[i] += u_right*((1.0)**k)*wenoCoeffs[k]
            f_left[i] += u_left*((-1.0)**k)*wenoCoeffs[k]

        f_right[i] *= right
        f_left[i] *= left

    residual = np.linspace(0,0,numCells+1)

    for i in range(len(f)):
        # integrate the fluxes
        flux = f_right[i]*dt/dx
        residual[i+1] += flux
        residual[i] -= flux

    for i in range(len(f)):
        f_out[i] += residual[i]
        
    residual = np.linspace(0,0,numCells+1)
        
    for i in range(len(f)):
        # integrate the fluxes
        flux = - f_left[i]*dt/dx
        residual[i] += flux
        residual[i+1] -= flux

    for i in range(len(f)):
        f_out[i] += residual[i+1]
    
    return f_out

