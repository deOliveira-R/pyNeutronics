import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from domain import Domain
import diffusion_fdm as dif

sns.set()


def dif_func(r):
    value = 5.0 * (r <= 5) \
          + 1.0 * (5 < r)
    return value


def sig_a_func(r):
    value = 1.0 * (r <= 5) \
          + 0.1 * (5 < r)
    return value


def nu_fission_func(r):
    value = 0.4 * (r <= 5) \
          + 0.0 * (5 < r)
    return value


def source_func(r):
    value = 0 * (r <= 5) \
          + 1.0 * (5 < r)
    return value


R = 10
I = 30

BC = (1, 0, 0)  # flux 0 at boundary

slab = Domain.uniform(R, I, 'slb')
slab_diffusion = dif.DiffusionSource.from_position_function(slab, dif_func, sig_a_func, nu_fission_func, source_func)
phi_slb = slab_diffusion.solve(BC)


def slab_ana(r):
    value = (r <= 5) * 0.181299 * np.exp(-0.34641 * r) * (np.exp(0.69282 * r) + 1) \
          + (5 < r) * (34.5868 * np.sinh(0.316228 * r) - 35.3082 * np.cosh(0.316228 * r) + 10)
    return value


phi_ana = slab_ana(slab.centers)

cylinder = Domain.uniform(R, I, 'cyl')
cylinder_diffusion = dif.DiffusionSource.from_position_function(cylinder, dif_func, sig_a_func, nu_fission_func, source_func)
phi_cyl = cylinder_diffusion.solve(BC)

sphere = Domain.uniform(R, I, 'sph')
sphere_diffusion = dif.DiffusionSource.from_position_function(sphere, dif_func, sig_a_func, nu_fission_func, source_func)
phi_sph = sphere_diffusion.solve(BC)

plt.figure()
ax = plt.plot(slab.centers, phi_slb)
ax_ana = plt.scatter(slab.centers, phi_ana)
ax2 = plt.plot(cylinder.centers, phi_cyl)
ax3 = plt.plot(sphere.centers, phi_sph)
plt.show()

# Reed's problem


def dif_func(r):
    value = 1.0/3.0 * (5 < r) \
          + 1.0/3.0/0.001 * ((3 < r) * (r <= 5)) \
          + 1.0/3.0/5.0 * ((2 < r) * (r <= 3)) \
          + 1.0/3.0/50.0 * (r <= 2)
    return value


def sig_a_func(r):
    value = 0\
          + (0.1 * (5 < r)
          + 5.0 * ((2 < r) * (r < 3))
          + 50.0 * (r <= 2))
    return value


def nu_fission_func(r):
    return 0*r


def source_func(r):
    value = 0\
          + 1.0 * ((5 < r) * (r < 7))\
          + 50.0 * (r <= 2)
    return value


R = 9
I = 50
BC = (0.25, 0.5 * dif_func(R), 0)  # vacuum BC

slab = Domain.uniform(R, I, 'slb')
slab_diffusion = dif.DiffusionSource.from_position_function(slab, dif_func, sig_a_func, nu_fission_func, source_func)
phi_slb = slab_diffusion.solve(BC)

cylinder = Domain.uniform(R, I, 'cyl')
cylinder_diffusion = dif.DiffusionSource.from_position_function(cylinder, dif_func, sig_a_func, nu_fission_func, source_func)
phi_cyl = cylinder_diffusion.solve(BC)

sphere = Domain.uniform(R, I, 'sph')
sphere_diffusion = dif.DiffusionSource.from_position_function(sphere, dif_func, sig_a_func, nu_fission_func, source_func)
phi_sph = sphere_diffusion.solve(BC)

plt.figure()
ax = plt.plot(slab.centers, phi_slb)
ax2 = plt.plot(cylinder.centers, phi_cyl)
ax3 = plt.plot(sphere.centers, phi_sph)
plt.show()
