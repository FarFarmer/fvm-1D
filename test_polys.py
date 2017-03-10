import numpy as np
import matplotlib.pyplot as plt
import fvm

polyOrder = 2
spaceOrder = polyOrder+1
res = 10
f = np.linspace(0,0,res)
f[2] = 1.0
f[3] = 0.75
f[4] = 0.85

# # this shows problems
# f[0] = 0
# f[1] = 0
# f[2] = 8.52760128e-02
# f[3] = 5.11621972e-01                                                
# f[4] = 3.96826532e-01

f[0] = 5.2 #/5.5;
f[1] = 5.5 #/5.5;
f[2] = 5.0 #/5.5;
f[3] = 1.3 #/5.5;
f[4] = 1.1 #/5.5;
f[5] = 0.8 #/5.5;

  
u = np.linspace(1,1,res)
dx = 2
k = 2 # central cell

L = fvm.generateStencilMat(polyOrder)
Linv = np.linalg.inv(L)
polyCoeffs = fvm.calcPolyCoeffs(polyOrder, f, k, Linv)
S = fvm.generateOscillMat(polyOrder, dx)
omegas = fvm.calcOmega(polyCoeffs, S, polyOrder)
wenoCoeffs = np.transpose(polyCoeffs).dot(omegas)

# print("polyCoeffs:");
# print(polyCoeffs);
# print("wenoCoeffs:")
# print(wenoCoeffs)
# print("omegas:")
# print(omegas)

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
    
    for p in range(spaceOrder):
        print("w = " + str(polyCoeffs[polyIdx][p]) +" offset = " + str(x_k) + " scale = " + str(fvm.cellWidth/dx) + " p = " + str(p))

    print("")

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

print("weno coeffs")
for p in range(spaceOrder):
    print("w = " + str(wenoCoeffs[p]) +" offset = " + str(x_k) + " scale = " + str(fvm.cellWidth/dx) + " p = " + str(p))    
        
plt.plot(x, y, '--', linewidth=2)

y_right_k_face = 0
for p in range(spaceOrder):
    y_right_k_face += wenoCoeffs[p]*((dx/2)*fvm.cellWidth/dx)**p


plt.plot(np.linspace(-0.5*dx,(res-0.5)*dx,2),y_right_k_face*np.linspace(1,1,2),
         linestyle='--', color='black')   
    
ind = np.arange(-0.5*dx,(res-0.5)*dx,dx)
plt.bar(ind, f, dx, fill=False)
plt.xlabel('x')
plt.ylabel('f')
plt.title('reconstruction')
plt.grid(True)

plt.show()
