import numpy as np
import matplotlib.pyplot as plt
from shocktubecalc import sod

N = 1000

positions, regions, Q_exact = sod.solve(left_state=(1., 1., 0.), right_state=(0.1, 0.125, 0.), geometry=(0., 1., 0.5), t=0.15, gamma=1.4, npts=N)

print(positions)

gamma = 1.4

Q0_euler = np.genfromtxt('sod_shock_initial_conditions_euler.csv', delimiter = ',')
Qf_e_mac = np.genfromtxt('sod_shock_final_euler_maccormack.csv', delimiter = ',')
Qf_e_musc = np.genfromtxt('sod_shock_final_euler_muscl.csv', delimiter = ',')

rho0 = Q0_euler[0]
u0 = Q0_euler[1]/Q0_euler[0]
e0 = Q0_euler[2]

rhof_mac = Qf_e_mac[0]
uf_mac = Qf_e_mac[1]/Qf_e_mac[0]
ef_mac = Qf_e_mac[2]

rhof_mus = Qf_e_musc[0]
uf_mus = Qf_e_musc[1]/Qf_e_musc[0]
ef_mus = Qf_e_musc[2]

rhoy = Q_exact['rho']
uy = Q_exact['u']
ey = (1/(gamma - 1))*Q_exact['p'] + 0.5*rhoy*uy**2

#Q0_mhd = np.genfromtxt('sod_shock_initial_conditions_mhd.csv', delimiter = ',')
#Qf_m_mac = np.genfromtxt('sod_shock_final_mhd_maccormack.csv', delimiter = ',')
#Qf_m_musc = np.genfromtxt('sod_shock_final_mhd_muscl.csv', delimiter = ',')

error_rho_mac = np.sqrt((1/N)*np.sum((rhoy - rhof_mac)**2))
error_rho_mus = np.sqrt((1/N)*np.sum((rhoy - rhof_mus)**2))
error_u_mac = np.sqrt((1/N)*np.sum((uy - uf_mac)**2))
error_u_mus = np.sqrt((1/N)*np.sum((uy - uf_mus)**2))
error_e_mac = np.sqrt((1/N)*np.sum((ey - ef_mac)**2))
error_e_mus = np.sqrt((1/N)*np.sum((ey - ef_mus)**2))

print(error_rho_mac)
print(error_rho_mus)
print(error_u_mac)
print(error_u_mus)
print(error_e_mac)
print(error_e_mus)

x = np.linspace(0, 1, num=len(Q0_euler[0]))

plt.figure()
plt.plot(x, rhof_mac, label="MacCormack")
plt.plot(x, rhof_mus, label="MUSCL")
plt.plot(x, rhoy, label="Exact")
plt.xlabel("$x$")
plt.ylabel(r"$\rho$")
plt.legend()
plt.show()

plt.figure()
plt.plot(x, uf_mac, label="MacCormack")
plt.plot(x, uf_mus, label="MUSCL")
plt.plot(x, uy, label="Exact")
plt.xlabel("$x$")
plt.ylabel("$u$")
plt.legend()
plt.show()

plt.figure()
plt.plot(x, ef_mac, label="MacCormack")
plt.plot(x, ef_mus, label="MUSCL")
plt.plot(x, ey, label="Exact")
plt.xlabel("$x$")
plt.ylabel("$P$")
plt.legend()
plt.show()
