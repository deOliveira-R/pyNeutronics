import numpy as np
from scipy.special import j0
import matplotlib.pyplot as plt
import seaborn as sns
import diffusion_fdm as dif

sns.set()


def dif_func(r):
    return 4.0*(r > 0)  # np.full_like(r, 4)


def sig_a_func(r):
    return np.full_like(r, 0.18)


def nu_fission_func(r):
    return np.full_like(r, 0.2)


# Domain definition
R = 100
I = 40

dif_coe = dif_func(R)

BC = (0.25, 0.5 * dif_func(R))
epsilon = 1.0e-10

slab = dif.Domain.uniform(R, I, 'slab')
slab_diffusion = dif.DiffusionEigenvalue1E.from_position_function(slab, dif_func, sig_a_func, nu_fission_func)
k, phi_slb = slab_diffusion.solve(BC, epsilon=epsilon)

Re = R + 0.7104*3*dif_coe
Bg = np.pi / (2*Re)
slb_ana = np.cos(Bg * slab.centers)
slb_ana /= np.linalg.norm(slb_ana)

cylinder = dif.Domain.uniform(R, I, 'cyl')
cylinder_diffusion = dif.DiffusionEigenvalue1E.from_position_function(cylinder, dif_func, sig_a_func, nu_fission_func)
kc, phi_cyl = cylinder_diffusion.solve(BC, epsilon=epsilon)

Bg = 2.405 / Re
cyl_ana = j0(Bg * cylinder.centers)
cyl_ana /= np.linalg.norm(cyl_ana)

sphere = dif.Domain.uniform(R, I, 'sph')
sphere_diffusion = dif.DiffusionEigenvalue1E.from_position_function(sphere, dif_func, sig_a_func, nu_fission_func)
ks, phi_sph = sphere_diffusion.solve(BC, epsilon=epsilon)

Bg = np.pi / Re
sph_ana = np.sin(Bg * sphere.centers)/sphere.centers
sph_ana /= np.linalg.norm(sph_ana)

plt.figure()
ax = plt.plot(slab.centers, phi_slb)
# axa = plt.scatter(slab.centers, slb_ana)
ax2 = plt.plot(cylinder.centers, phi_cyl)
# axa = plt.scatter(cylinder.centers, cyl_ana)
ax3 = plt.plot(sphere.centers, phi_sph)
# axa  = plt.scatter(sphere.centers, sph_ana)
plt.show()
