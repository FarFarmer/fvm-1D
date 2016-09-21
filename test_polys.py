import numpy as np
import matplotlib.pyplot as plt
import fvm

polyOrder = 2
res = 6
f = np.linspace(0,0,res)
f[2] = 1.0
f[3] = 1.0

u = np.linspace(1,1,res)

dx = 1
dt = 0.125

Linv = np.linalg.inv(generateStencilMat(polyOrder))
S = generateOscillMat(polyOrder, dx)

polyCoeffs = np.ndarray(shape=(polyOrder+1,polyOrder+1))
for j in range(polyOrder+1):
    
    # ids is used again to avoid if-else blocks regarding values on
    # the boundary; it clamps the indices to [0:numCells-1]
    ids = list(range(i-polyOrder+j,i+j+1))
    for k in range(len(ids)):
        ids[k] = max(0, min(ids[k], numCells-1))
        
        # scale = 2 because in the reference frame, dx = 2
        scale = 2.0
        phi = np.array([scale*f[x] for x in ids])
        polyCoeffs[j] = Linv[j].dot(phi)
