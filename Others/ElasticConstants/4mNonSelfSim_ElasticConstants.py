# Elastic constants used for the analysis

import numpy as np

# Metagabbro

# 1. From waveform inversion of coupling effect
cp = 6200 #[m/s]
cs = 3600
rho = 2980

nu = 0.5 * ((cp/cs)**2-2)/((cp/cs)**2-1)

mu = rho * cs**2
E = 2* mu * (1+nu)

print("Waveform inv: from cp, cs and rho.")
print(f"E={E/1e9:.5f} GPa")
print(f"mu={mu/1e9:.5f} GPa")
print(f"nu={nu:.5f}")
print(f"vp={cp} m/s")
print(f"vs={cs} m/s")
print(f"rho={rho} kg/m^3")


# 2. used for the dynamic rupture simulation

# Elastic constant
print("Dynamic rupture model: from E, nu, rho")
E_dyn = 96e9
rho_dyn = 2980
nu_dyn = 0.246 # metagabbro
mu_dyn = E_dyn/(2*(1+nu_dyn))
cs_dyn = np.sqrt(mu_dyn/rho_dyn)
cp_dyn = cs_dyn * np.sqrt(2*(nu-1)/(2*nu-1))

print(f"E={E_dyn/1e9:.5f} GPa")
print(f"mu={mu_dyn/1e9:.5f} GPa")
print(f"nu={nu_dyn:.5f}")
print(f"vp={cp_dyn:.3f} m/s")
print(f"vs={cs_dyn:.3f} m/s")
print(f"rho={rho_dyn} kg/m^3")


"""
Output:
$python 4mNonSelfSim_ElasticConstants.py

Waveform inv: from cp, cs and rho.
E=96.21854 GPa
mu=38.62080 GPa
nu=0.24568
vp=6200 m/s
vs=3600 m/s
rho=2980 kg/m^3
Dynamic rupture model: from E, nu, rho
E=96.00000 GPa
mu=38.52327 GPa
nu=0.24600
vp=6192.167 m/s
vs=3595.452 m/s
rho=2980 kg/m^3
"""

