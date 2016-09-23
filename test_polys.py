import numpy as np
import matplotlib.pyplot as plt
import fvm

polyOrder = 2
spaceOrder = polyOrder+1
res = 6
f = np.linspace(0,0,res)
f[2] = 1.0
f[3] = 1.0

u = np.linspace(1,1,res)
dx = 1
org = -0.5
k = 3 # central cell

L = fvm.generateStencilMat(polyOrder)
Linv = np.linalg.inv(L)
polyCoeffs = fvm.calcPolyCoeffs(polyOrder, f, k, Linv)
S = fvm.generateOscillMat(polyOrder, dx)
omegas = fvm.calcOmega(polyCoeffs, S, polyOrder)
wenoCoeffs = np.transpose(polyCoeffs).dot(omegas)

polyRes = 101

polyIdxExtent = int(polyOrder/2+0.5)
polyExtent = (polyIdxExtent+0.5)*dx

polyIdx = 0
for i in range(k-polyIdxExtent,k+polyIdxExtent+1):

    # p-space values
    x0 = -polyExtent+i*dx # start x-coordinate of the polynomial
    x1 =  polyExtent+i*dx # end x-coordinate of the polynomial
    x_k = k*dx
    
    x = np.linspace(x0,x1,polyRes)
    y = np.linspace(0,0,polyRes)

    for j in range(polyRes):
        for p in range(spaceOrder):
            y[j] += polyCoeffs[polyIdx][p]*((x[j]-x_k)*fvm.cellWidth/dx)**p
    
    plt.plot(x, y)
    polyIdx += 1

    
# p-space values
x0 = -polyExtent+k*dx # start x-coordinate of the polynomial
x1 =  polyExtent+k*dx # end x-coordinate of the polynomial
x_k = k*dx

x = np.linspace(x0,x1,polyRes)
y = np.linspace(0,0,polyRes)

for j in range(polyRes):
    for p in range(spaceOrder):
        y[j] += wenoCoeffs[p]*((x[j]-x_k)*fvm.cellWidth/dx)**p
        
plt.plot(x, y, '--', linewidth=2)


    

ind = np.arange(-0.5*dx,(res-0.5)*dx,dx)
plt.bar(ind, f, dx, color='r')
plt.xlabel('x')
plt.ylabel('f')
plt.title('reconstruction')
plt.grid(True)
    
plt.show()
