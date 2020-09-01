import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import diffusion_fdm as dif

sns.set()


def dif_func(r):
    return 5.0*(r <= 5) + 1.0*(r > 5)


def sig_a_func(r):
    return 0.5*(r <= 5) + 0.01*(r > 5)


def nu_fission_func(r):
    return 0.7*(r <= 5) + 0.0*(r > 5)


R = 10
I = 10

BC = (1, 0)

slab = dif.Domain.uniform(R, I, 'slab')
cylinder = dif.Domain.uniform(R, I, 'cylindrical')
sphere = dif.Domain.uniform(R, I, 'spherical')

slab_diffusion = dif.DiffusionEigenvalue1E.from_position_function(slab, dif_func, sig_a_func, nu_fission_func)
cylinder_diffusion = dif.DiffusionEigenvalue1E.from_position_function(cylinder, dif_func, sig_a_func, nu_fission_func)
sphere_diffusion = dif.DiffusionEigenvalue1E.from_position_function(sphere, dif_func, sig_a_func, nu_fission_func)

k,  phi_slb = slab_diffusion.solve(BC)
kc, phi_cyl = cylinder_diffusion.solve(BC)
ks, phi_sph = sphere_diffusion.solve(BC)

plt.figure()
ax = plt.plot(slab.centers, phi_slb)
ax2 = plt.plot(cylinder.centers, phi_cyl)
ax3 = plt.plot(sphere.centers, phi_sph)
plt.show()


def dif_func(r):
    return 1.16356040*(r <= 75)\
           + 2.04235228*(r > 75)


def sig_a_func(r):
    return 0.0354966*(r <= 75)\
           + 0.01666126*(r > 75)


def nu_fission_func(r):
    return 0.036144338*(r <= 75)\
           + 0.0*(r > 75)


# d_extr = 0.7104*dif_func(100)
#
# R = 275
# I = 500
#
# L = lambda D, A: np.sqrt(D/A)
# Lf = L(dif_func(5), sig_a_func(5))
# Lr = L(dif_func(100), sig_a_func(100))
#
# C = dif_func(100)/Lr/np.tanh(200+d_extr)
#
# Bg = 0.019805
# keff = 1.005859
# Bm = np.sqrt((nu_fission_func(5)/keff - sig_a_func(5)) / dif_func(5))
#
# A = np.cos(Bm * 75)/np.sinh((200 + d_extr)/Lr)
# slab_fuel_ana = lambda Bm, x: np.cos(Bm * x)
# slab_refl_ana = lambda A, Re, Lr, x: A * np.sinh((Re - x)/Lr)
#
#
# slab = dif.Domain.uniform(R, I, 'slab')
#
# slab_diffusion = dif.DiffusionEigenvalue1E.from_position_function(slab, dif_func, sig_a_func, nu_fission_func)
# # kn, phi_num, centers = slab_diffusion.solve(BC)
#
# slab_ana = lambda r: slab_fuel_ana(Bm, r) * (r <= 75) \
#                    + slab_refl_ana(A, R + d_extr, Lr, r) * (r > 75)
#
# slab_res = slab_ana(slab.centers)
# slab_res /= np.linalg.norm(slab_res)  # make norm(x)==1

# plt.figure()
# axa = plt.plot(slab.centers, phi_num)
# axb = plt.plot(slab.centers, slab_res)
# axb = plt.plot(slab.centers, slab_refl_ana_res)
# plt.show()