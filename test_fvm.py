import numpy as np
import matplotlib.pyplot as plt
import fvm

polyOrder = 2

res = 100 #input

f = np.linspace(0,0,res)
f[2] = 1.0

u = np.linspace(1,1,res)

dx = 0.5
dt = 0.25

N = res
ind = np.arange(-0.5,N-0.5,1)  # the x locations for the groups
width = 1       # the width of the bars

fig, ax = plt.subplots()

plt.ion()
plt.show()

def func(coeffs,x):
    scale = 2.0
    # print(coeffs)
    y = np.linspace(0,0,len(x))
    for i in range(len(x)):
        y[i] = coeffs[2]*(x[i]*scale)**2 + coeffs[1]*(x[i]*scale) + coeffs[0]
    return y

while True:
    
    ax.clear()    
    ax.set_ylabel('f(x)')
    ax.set_xlabel('x')
    ax.set_title('Dye advection')
    ax.set_ylim(top=1.0,bottom=-0.2)
    ax.set_xlim(left=0,right=N)
    # ax.set_xticks(ind + width/2.0)
    # ax.set_xticklabels([str(idx) for idx in range(N)])

    f_new = fvm.advectWENO(f,u,dt,dx,polyOrder)
    # f_new = fvm.advect(f,u,dt,dx)
    rects = ax.bar(ind, f_new, width, color='r')

    # # test polys
    # x = np.arange(-N,N,0.1)
    # y = func(wenoCoeffs,x)
    # plt.plot(x,y)
    # #~test polys
    
    plt.draw()
    
    f,f_new = f_new,f
    
    input()
