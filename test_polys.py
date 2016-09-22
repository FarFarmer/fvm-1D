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
dx = 2
k = 3

L = fvm.generateStencilMat(polyOrder)
Linv = np.linalg.inv(L)
S = fvm.generateOscillMat(polyOrder, dx)

polyCoeffs = fvm.calcPolyCoeffs(polyOrder, f, k, Linv)

#for each polynomial
for i in range(-int(polyOrder/2),int(polyOrder/2)+1)
    idx = 

# omegas = fvm.calcOmega(polyCoeffs, S, polyOrder)
# wenoCoeffs = np.transpose(polyCoeffs).dot(omegas)


print(wenoCoeffs)


x = np.linspace(-1,1,100)
y = np.linspace(0,1,100)

for i in range(len(x)):
    y[i] = 0
    for j in range(polyOrder):
        y[i] += wenoCoeffs[j]*x[i]**j

plt.plot(x, y)

ind = np.arange(-0.5,res-0.5,1)
width = 1
plt.bar(ind, f, width, color='r')
#plt.xlim(left=-0.5,right=res-0.5)
plt.xlim(left=-1.5,right=res-0.5)
plt.xlabel('x')
plt.ylabel('f')
plt.title('reconstruction')
plt.grid(True)
plt.show()



