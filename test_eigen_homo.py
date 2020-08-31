import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import diffusion_fdm as dif

sns.set()


def dif_func(r):
    return np.full_like(r, 4)


def sig_a_func(r):
    return np.full_like(r, 0.18)


def nu_fission_func(r):
    return np.full_like(r, 0.2)


# Domain definition
R = 100
I = 10

BC = (1, 0)
epsilon = 1.0e-10

slab = dif.Domain.uniform(R, I, 'slab')
slab_diffusion = dif.DiffusionEigenvalue1E.from_position_function(slab, dif_func, sig_a_func, nu_fission_func)
k_new, phi_slb, centers = slab_diffusion.solve(BC, epsilon=epsilon)
# phi_slb_ana = A * np.cos(np.pi*centers/R)

cylinder = dif.Domain.uniform(R, I, 'cyl')
cylinder_diffusion = dif.DiffusionEigenvalue1E.from_position_function(cylinder, dif_func, sig_a_func, nu_fission_func)
kc, phi_cyl, centers = cylinder_diffusion.solve(BC, epsilon=epsilon)

sphere = dif.Domain.uniform(R, I, 'sph')
sphere_diffusion = dif.DiffusionEigenvalue1E.from_position_function(sphere, dif_func, sig_a_func, nu_fission_func)
ks, phi_sph, centers = sphere_diffusion.solve(BC, epsilon=epsilon)

plt.figure()
ax = plt.plot(centers, phi_slb)
ax2 = plt.plot(centers, phi_cyl)
ax3 = plt.plot(centers, phi_sph)
plt.show()

# solution in cylindrical geometry with 20 cells
# A_cyl_orig, B_cyl_orig, k_c_orig, phi_cyl_orig, centers = dif.DiffusionEigenvalueOriginal(R, I, dif_func, sig_a_func, nu_fission_func, BC, 'cyl', epsilon=epsilon)
# phi_cyl_ana = C * j0(2.405*centers/R)

# solution in spherical geometry with 20 cells
# A_sph_orig, B_sph_orig, k_s_orig, phi_sph_orig, centers = dif.DiffusionEigenvalueOriginal(R, I, dif_func, sig_a_func, nu_fission_func, BC, 'sph', epsilon=epsilon)
# phi_sph_ana = A/centers * np.sin(np.pi*centers/R)
