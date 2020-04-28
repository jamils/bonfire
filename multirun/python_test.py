import numpy as np
import matplotlib.pyplot as plt
from shocktubecalc import sod

N1 = 1000
N2 = 100
N3 = 10

positions1, regions1, Q1 = sod.solve(left_state=(1., 1., 0.), right_state=(0.1, 0.125, 0.), geometry=(0., 1., 0.5), t=0.15, gamma=1.4, npts=N1)
positions2, regions2, Q2 = sod.solve(left_state=(1., 1., 0.), right_state=(0.1, 0.125, 0.), geometry=(0., 1., 0.5), t=0.15, gamma=1.4, npts=N2)
positions3, regions3, Q3 = sod.solve(left_state=(1., 1., 0.), right_state=(0.1, 0.125, 0.), geometry=(0., 1., 0.5), t=0.15, gamma=1.4, npts=N3)

print(positions1)

rho1 = Q1['rho']
rho2 = Q2['rho']
rho3 = Q3['rho']

x1 = np.linspace(0, 1, num=len(rho1))
x2 = np.linspace(0, 1, num=len(rho2))
x3 = np.linspace(0, 1, num=len(rho3))

plt.figure()
plt.scatter(x1, rho1, label="N = 1000")
plt.scatter(x2, rho2, label="N = 100")
plt.scatter(x3, rho3, label="N = 10")
plt.xlabel("$x$")
plt.ylabel(r"$\rho$")
plt.legend()
plt.show()
