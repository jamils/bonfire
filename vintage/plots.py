import numpy as np
import matplotlib.pyplot as plt

Q0 = np.genfromtxt('sod_shock_initial_conditions_mhd.csv', delimiter = ',')
Qf = np.genfromtxt('sod_shock_final_mhd_maccormack.csv', delimiter = ',')

N = len(Q0[0,:])

rho0 = Q0[0]
u0 = Q0[1]/Q0[0]
E0 = Q0[2]
x = np.linspace(0, 1, num=N)

rhof = Qf[0]
uf = Qf[1]/Qf[0]
Ef = Qf[2]
x = np.linspace(0, 1, num=N)

plt.plot(x, rho0, label="$\rho_0$")
plt.plot(x, rhof, label="$\rho_f$")
plt.legend()
plt.show()
