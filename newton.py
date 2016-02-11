# Plotkin, Benjamin
# 18th April, 2015
# Newtonian Fractal Generator

import math
#import cmath as cm
import numpy as np
import matplotlib.pyplot as plt

# finds nth roots of unity
def unityRoot(n):
    #roots = np.zeros([n], dtype=complex)
    roots = np.zeros(n,complex)
    for k in range(n):
        # prof's way (no cmath)
        theta = 2*math.pi*k/n
        roots[k] = math.cos(theta)+math.sin(theta)*1j
        # my way (uses cmath)
        #roots[k] = cm.exp(cm.sqrt(-1)*2*cm.pi*k/n)
    return roots

# generates the newton fractal for a given polynomial when called from the
# interpreter command line -- good starting values are as follows:
# K = 30, pixel = 400, tol = .001, n = 2, 3, or 4
def nfractal(K, pixel, tol, n):
    r = 0.75; # ratio of number points along imaginary axis to number of points
              # along real axis
    a = 1     # 0.75 - 0.4j #1.0 #0.4;

    x_left = -1.5;
    x_right = 1.5;
    y_bottom = -1.5;
    y_top = 1.5;

    # x / real coordinates
    x = np.linspace(x_left, x_right, pixel)
    # y / imaginary coordinates
    y = np.linspace(y_bottom, y_top ,round(pixel*r))

    # meshgrid
    [Re, Im] = np.meshgrid(x,y)
    C = np.zeros([round(r*pixel), pixel], complex)
    # matrix holding x and y coordinates or points on complex plane
    C = Re + Im*1.0j
    # this matrix keeps track of speed of convergence to a root
    # if B(i,j) is large, iteration diverges...
    B = np.zeros([round(r*pixel), pixel], float)
    # matrix of ones
    Id = np.ones(B.shape);

    nRoots=unityRoot(n)

    # initial guess for Newton's iteration
    Cn = C

    # iterate the initial guess k times; the iteration is applied to the whole
    # matrix using array element-by-element operations
    for k in range(1,K):
        Cn -= a*(Cn**n - Id) / (n*Cn**(n-1))
        for j in xrange(0,n):
            # increment B by one if the corresponding entry in matrix Cn
            # remains close to a root of the polynomial after the kth iteration
            B += (np.abs(Cn-Id*nRoots[j] < tol))

    # plot fractal
    plt.pcolormesh(x, y, B, cmap = 'bone')
    plt.colorbar()
    plt.show()
