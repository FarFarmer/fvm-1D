import numpy as np
import matplotlib.pyplot as plt
import fvm
import time

polyOrder = 4

res = 128 #input

u = np.linspace(0,0,res)
u[4] = 1.0

v = np.linspace(1,1,res)

dx = 0.25
dt = 0.125

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
    ax.set_ylabel('u(x)')
    ax.set_xlabel('x')
    ax.set_title('Dye advection')
    ax.set_ylim(top=1.0,bottom=-0.2)
    ax.set_xlim(left=0,right=N)
    # ax.set_xticks(ind + width/2.0)
    # ax.set_xticklabels([str(idx) for idx in range(N)])

    
    u_new = fvm.advectWENO(u,v,dt,dx,polyOrder)
    # u_new = fvm.advect(u,v,dt,dx)
    rects = ax.bar(ind, u_new, width, color='r')

    # # test polys
    # x = np.arange(-N,N,0.1)
    # y = func(wenoCoeffs,x)
    # plt.plot(x,y)
    # #~test polys
    
    plt.draw()
    
    u,u_new = u_new,u
    
    input()
