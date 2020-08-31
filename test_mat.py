import xs_library_1e as lib
from material_data import Material, NuclearMaterial
import diffusion_fdm as dif

import matplotlib.pyplot as plt
import seaborn as sns

sns.set()


def mapping_homo(r):
    if r <= 100:
        return homogeneous_fuel


def mapping_hetero(r):
    if r <= 80:
        return homogeneous_fuel
    elif 100 < r:
        return moderator


R = 100
I = 10

BC = [1, 0]
epsilon = 1.0e-10

slab = dif.Domain.uniform(R, I, 'slab')

#### MAKE INDIVIDUAL MATERIALS

fuel_enrichment = 0.035  # 3.5% U235  critical at 1.65%
moderator_to_fuel_ratio = 1.668  # original 1.668   critical at 16.35

fuel = Material(-10.4, {'U235': fuel_enrichment,
                        'U238': 1 - fuel_enrichment,
                        'O16': 2})

moderator = Material(-0.75, {'H1': 2, 'O16': 1})

homogeneous_fuel = Material.mix(moderator, fuel, moderator_to_fuel_ratio)

### GIVE XS TO MATERIAL (Make it into a material with nuclear properties!)

homogeneous_fuel = NuclearMaterial(lib, homogeneous_fuel)

materials = list(map(mapping_homo, slab.centers))

print(materials)

slab_diffusion = dif.DiffusionEigenvalue1E.from_materials(slab, materials)
k, phi, centers = slab_diffusion.solve(BC, epsilon=epsilon)
# phi_slb_ana = A * np.cos(np.pi*centers/R)

plt.figure()
ax = plt.plot(centers, phi)
plt.show()
