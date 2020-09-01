import numpy as np
from scipy.special import k0
import matplotlib.pyplot as plt
import seaborn as sns

from domain import Domain
import diffusion_fdm as dif

sns.set()


def dif_func(r):
    return np.full_like(r, 0.04)


def sig_a_func(r):
    return np.full_like(r, 1)


def nu_fission_func(r):
    return np.full_like(r, 0)


def source_func(r):
    return np.full_like(r, 1)


print("For this problem the diffusion length is", np.sqrt(dif_func(1)/sig_a_func(1)))
inf_med = source_func(1)/(sig_a_func(1) - nu_fission_func(1))
print("The infinite medium solution is", inf_med)

R = 10
I = 10

BC = (0, 1, 0)  # reflective BC

slab = Domain.uniform(R, I, 'slb')
cylinder = Domain.uniform(R, I, 'cyl')
sphere = Domain.uniform(R, I, 'sph')

slab_diffusion = dif.DiffusionSource.from_position_function(slab, dif_func, sig_a_func, nu_fission_func, source_func)
cylinder_diffusion = dif.DiffusionSource.from_position_function(cylinder, dif_func, sig_a_func, nu_fission_func, source_func)
sphere_diffusion = dif.DiffusionSource.from_position_function(sphere, dif_func, sig_a_func, nu_fission_func, source_func)

phi_slb = slab_diffusion.solve(BC)
phi_cyl = cylinder_diffusion.solve(BC)
phi_sph = sphere_diffusion.solve(BC)

plt.figure()
ax = plt.plot(slab.centers, phi_slb)
ax2 = plt.plot(cylinder.centers, phi_cyl)
ax3 = plt.plot(sphere.centers, phi_sph)
plt.show()


def dif_func(r):
    return 0.4 * (0.0 < r)


def sig_a_func(r):
    return 0.1 * (0.0 < r)


def nu_fission_func(r):
    return np.full_like(r, 0)


# def source_func(r):
#     return np.full_like(r, 1)

def source_func(r):
    value = 2.0 * (r <= 0.02) \
          + 0 * (0.02 < r)
    return value

R = 10
I = 500

source_pos = R/I
const_pos = R - 1

S0 = source_func(0.01)* 0.02 * 2

dif_coe = dif_func(5)
sig_a_coe = sig_a_func(5)
L = np.sqrt(dif_coe/sig_a_coe)

C_slb = S0*L / (2*dif_coe)
C_cyl = S0/(2*dif_coe*np.pi)
C_sph = S0/(4*dif_coe*np.pi)

BC = (0.25, 0.5 * dif_func(R), 0)  # vacuum BC

slab = Domain.uniform(R, I, 'slb')
cylinder = Domain.uniform(R, I, 'cyl')
sphere = Domain.uniform(R, I, 'sph')

slab_diffusion = dif.DiffusionSource.from_position_function(slab, dif_func, sig_a_func, nu_fission_func, source_func)
cylinder_diffusion = dif.DiffusionSource.from_position_function(cylinder, dif_func, sig_a_func, nu_fission_func, source_func)
sphere_diffusion = dif.DiffusionSource.from_position_function(sphere, dif_func, sig_a_func, nu_fission_func, source_func)

phi_slb = slab_diffusion.solve(BC)
slab_ana = lambda C, L, r: C * np.exp(-r/L)
slab_ana_res = slab_ana(C_slb, L, slab.centers)

phi_cyl = cylinder_diffusion.solve(BC)
cyl_ana = lambda C, L, r: C * k0(r/L)
cyl_ana_res = cyl_ana(C_cyl, L, cylinder.centers)

phi_sph = sphere_diffusion.solve(BC)
sph_ana = lambda C, L, r: C * np.exp(-r/L)/r
sph_ana_res = sph_ana(C_sph, L, sphere.centers)

plt.figure()
# ax = plt.plot(slab.centers, phi_slb)
# ax = plt.scatter(slab.centers, slab_ana_res)
# ax2 = plt.plot(cylinder.centers, phi_cyl)
# ax = plt.scatter(cylinder.centers, cyl_ana_res)
ax3 = plt.plot(sphere.centers, phi_sph)
ax = plt.scatter(sphere.centers, sph_ana_res)
plt.show()
