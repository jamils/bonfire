import numpy as np
import matplotlib.pyplot as plt

Q0 = np.genfromtxt('mhd_wave_initial_conditions_mhd.csv', delimiter = ',')
Qf = np.genfromtxt('sod_shock_final_mhd_muscl.csv', delimiter = ',')

N = len(Q0[0,:])

rho0 = Q0[0]
u0 = Q0[1]/Q0[0]
v0 = Q0[2]/Q0[0]
w0 = Q0[3]/Q0[0]
Bx0 = Q0[4]
By0 = Q0[5]
Bz0 = Q0[6]
E0 = Q0[7]

rhof = Qf[0]
uf = Qf[1]/Qf[0]
vf = Qf[2]/Qf[0]
wf = Qf[3]/Qf[0]
Bxf = Qf[4]
Byf = Qf[5]
Bzf = Qf[6]
Ef = Qf[7]

error = np.sqrt((1/N)*np.sum((Qf - Q0)**2))
print(error)

x = np.linspace(0, 1, num=N)

#plt.plot(x, rho0, label="$rho_0$")
#plt.plot(x, rhof, label="$rho_f$")
#plt.plot(x, u0, label="$u_0$")
#plt.plot(x, uf, label="$u_f$")
#plt.plot(x, Bzf/np.sqrt(rhof), label="$u_f$")
#plt.plot(x, v0, label="$v_0$")
#plt.plot(x, vf, label="$v_f$")
plt.plot(x, w0, label="$w_0$")
plt.plot(x, wf, label="$w_f$")
#plt.plot(x, Bx0, label="$Bx_0$")
#plt.plot(x, Bxf, label="$Bx_f$")
#plt.plot(x, By0, label="$By_0$")
#plt.plot(x, Byf, label="$By_f$")
#plt.plot(x, Bz0, label="$Bz_0$")
#plt.plot(x, Bzf, label="$Bz_f$")
#plt.plot(x, E0, label="$E_0$")
#plt.plot(x, Ef, label="$E_f$")
plt.legend()
plt.show()