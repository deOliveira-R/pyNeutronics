import numpy as np
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

phi_slb, rx = slab_diffusion.solve(BC)
phi_cyl, rc = cylinder_diffusion.solve(BC)
phi_sph, rs = sphere_diffusion.solve(BC)

plt.figure()
ax = plt.plot(rx, phi_slb)
ax2 = plt.plot(rc, phi_cyl)
ax3 = plt.plot(rs, phi_sph)
plt.show()


# def dif_func(r):
#     return np.full_like(r, 0.04)
#
#
# def sig_a_func(r):
#     return np.full_like(r, 1)
#
#
# def nu_fission_func(r):
#     return np.full_like(r, 0)
#
#
# def source_func(r):
#     return np.full_like(r, 1)
#
#
# print("For this problem the diffusion length is", np.sqrt(dif_func(1)/sig_a_func(1)))
# inf_med = source_func(1)/(sig_a_func(1) - nu_fission_func(1))
# print("The infinite medium solution is", inf_med)
#
# R = 10
# I = 10
#
# BC = (0, 1, 0)  # reflective BC
#
# slab = Domain.uniform(R, I, 'slb')
# cylinder = Domain.uniform(R, I, 'cyl')
# sphere = Domain.uniform(R, I, 'sph')
#
# slab_diffusion = dif.DiffusionSource.from_position_function(slab, dif_func, sig_a_func, nu_fission_func, source_func)
# cylinder_diffusion = dif.DiffusionSource.from_position_function(cylinder, dif_func, sig_a_func, nu_fission_func, source_func)
# sphere_diffusion = dif.DiffusionSource.from_position_function(sphere, dif_func, sig_a_func, nu_fission_func, source_func)
#
# phi_slb, rx = slab_diffusion.solve(BC)
# phi_cyl, rc = cylinder_diffusion.solve(BC)
# phi_sph, rs = sphere_diffusion.solve(BC)
#
# plt.figure()
# ax = plt.plot(rx, phi_slb)
# ax2 = plt.plot(rc, phi_cyl)
# ax3 = plt.plot(rs, phi_sph)
# plt.show()
