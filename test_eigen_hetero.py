import numpy as np
import scipy.sparse as sps
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


d_extr = 0.7104*dif_func(100)

R = 275
I = 500

L = lambda D, A: np.sqrt(D/A)
Lf = L(dif_func(5), sig_a_func(5))
Lr = L(dif_func(100), sig_a_func(100))

C = dif_func(100)/Lr/np.tanh(200+d_extr)

Bg = 0.019805
keff = 1.005859
Bm = np.sqrt((nu_fission_func(5)/keff - sig_a_func(5)) / dif_func(5))

A = np.cos(Bm * 75)/np.sinh((200 + d_extr)/Lr)
slab_fuel_ana = lambda Bm, x: np.cos(Bm * x)
slab_refl_ana = lambda A, Re, Lr, x: A * np.sinh((Re - x)/Lr)


slab = dif.Domain.uniform(R, I, 'slab')

slab_diffusion = dif.DiffusionEigenvalue1E.from_position_function(slab, dif_func, sig_a_func, nu_fission_func)
kn, phi_num = slab_diffusion.solve(BC)

slab_ana = lambda r: slab_fuel_ana(Bm, r) * (r <= 75) \
                   + slab_refl_ana(A, R + d_extr, Lr, r) * (r > 75)

slab_res = slab_ana(slab.centers)
slab_res /= np.linalg.norm(slab_res)  # make norm(x)==1

plt.figure()
axa = plt.plot(slab.centers, phi_num)
axb = plt.plot(slab.centers, slab_res)
plt.show()


Fuel = 258.6291
CVessel = 261.1291
Blanket = 325.5525
BVessel = 328.0525
Reflector = 428.0525


dif_fuel = np.array([3.65870, 2.56102])
dif_refl = np.array([1.45225, 1.00673])


def dif_func(r):
    return np.tensordot(dif_fuel, (r <= 328.0525), axes=0) \
         + np.tensordot(dif_refl, (r > BVessel), axes=0)


def rem_func(r):
    return np.tensordot(np.array([2.49234E-3, 4.28807E-3]), (r <= 328.0525), axes=0) \
         + np.tensordot(np.array([1.11615E-3, 2.39873E-4]), (r > BVessel), axes=0)


def scatmat_func(r):
    return np.tensordot(np.array([[1.10133E-1, 0.0],
                                  [1.01102E-3, 1.38648E-1]]), (r <= 328.0525), axes=0) \
         + np.tensordot(np.array([[2.67519E-1, 0.0],
                                  [9.9763E-4, 3.326E-1]]), (r > BVessel), axes=0)


def nu_fission_func(r):
    return np.tensordot(np.array([2.215E-3, 2.61801E-3]), (r <= 328.0525), axes=0) \
         + np.tensordot(np.array([0.0, 0.0]), (r > BVessel), axes=0)


def chi_func(r):
    return np.tensordot(np.array([9.98476E-1, 1.52399E-3]), (r <= 328.0525), axes=0) \
         + np.tensordot(np.array([0.0, 0.0]), (r > BVessel), axes=0)


I = 500
BC = (1, 0)

cyl = dif.Domain.uniform(Reflector, I, 'cylindrical')

# test_fac = np.array([2, 1, 3])
#
# test_int = dif_func(76)
# # print(test_int)
#
# test_vec = dif_func(np.array([74, 75, 76]))
# # print(test_vec)
# # print(test_vec[0])
#
# test_scat = scatmat_func(np.array([74, 75, 76]))
# print(test_scat)
# print(test_scat*test_fac)

cyl_diffusionMG = dif.DiffusionEigenvalueMG.from_position_function(cyl, dif_func, rem_func, scatmat_func, nu_fission_func, chi_func)
# A, B = slab_diffusionMG.assemble(BC)
# print(A[0][0].toarray())
#
# A_block = sps.bmat(A).toarray()
# B_block = sps.bmat(B).toarray()

kmg, phimg = cyl_diffusionMG.solve(BC)

plt.figure()
axmg1 = plt.plot(cyl.centers, phimg[0])
axmg2 = plt.plot(cyl.centers, phimg[1])
plt.show()
