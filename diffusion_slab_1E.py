#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 00:59:33 2019

@author: rodrigo
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import nuclear_data_library_1E as lib

from material_data import *

#flux_average = 1E14

#### MAKE INDIVIDUAL MATERIALS

fuel_enrichment = 0.035  # 3.5% U235  critical at 1.65%
moderator_to_fuel_ratio = 1.668  # original 1.668   critical at 16.35

fuel = Material(-10.4, {'U235': fuel_enrichment,
                        'U238': 1 - fuel_enrichment,
                        'O16': 2})

moderator = Material(-0.75, {'H1':2, 'O16':1})

homogeneous_fuel = Material.mix(moderator, fuel, moderator_to_fuel_ratio)
        
### GIVE XS TO MATERIAL (Make it into a material with nuclear properties!)

homogeneous_fuel = NuclearMaterial(lib, homogeneous_fuel)

### SPATIAL DISCRETIZATION

slab_size = 100  # cm
delta = 1  # size of interval for discretization

mesh_i = int(slab_size/delta)
mesh_p = mesh_i + 1

# use linspace for time. Allow space to be nonuniform instead.
x = np.linspace(-slab_size/2, slab_size/2, mesh_p)

# Calculate gain operator

op_gain = np.zeros((mesh_p, mesh_p))

for i in range(mesh_p):
    if i!=0 and i!=mesh_p-1: # if not in boundaries
        op_gain[i,i+1] = 1/8*delta*homogeneous_fuel.nu_fission()
        op_gain[i,i]   = 3/8*delta*(homogeneous_fuel.nu_fission()*2)
        op_gain[i,i-1] = 1/8*delta*homogeneous_fuel.nu_fission()
    else:
        if i==0:
            op_gain[i,i+1] = 1/8*delta*homogeneous_fuel.nu_fission()
            op_gain[i,i]   = 3/8*delta*homogeneous_fuel.nu_fission()
        if i==mesh_p-1:
            op_gain[i,i]   = 3/8*delta*homogeneous_fuel.nu_fission()
            op_gain[i,i-1] = 1/8*delta*homogeneous_fuel.nu_fission()


# Calculate diffusion coefficient and diffusion operator

op_diff = np.zeros((mesh_p, mesh_p))

for i in range(mesh_p):
    if i!=0 and i!=mesh_p-1: # if not in boundaries
        op_diff[i,i+1] = -homogeneous_fuel.diffusion()/delta
        op_diff[i,i]   = (homogeneous_fuel.diffusion()*2)/delta
        op_diff[i,i-1] = -homogeneous_fuel.diffusion()/delta
    else:
        if i==0:
            op_diff[i,i+1] = -homogeneous_fuel.diffusion()/delta
            op_diff[i,i]   = homogeneous_fuel.diffusion()/delta + 0.5
        if i==mesh_p-1:
            op_diff[i,i]   = homogeneous_fuel.diffusion()/delta + 0.5
            op_diff[i,i-1] = -homogeneous_fuel.diffusion()/delta
        
# Calculate absorption operator

op_abs = np.zeros((mesh_p, mesh_p))

for i in range(mesh_p):
    if i!=0 and i!=mesh_p-1: # if not in boundaries
        op_abs[i,i+1] = 1/8*delta*homogeneous_fuel.absorption()
        op_abs[i,i]   = 3/8*delta*(homogeneous_fuel.absorption()*2)
        op_abs[i,i-1] = 1/8*delta*homogeneous_fuel.absorption()
    else:
        if i==0:
            op_abs[i,i+1] = 1/8*delta*homogeneous_fuel.absorption()
            op_abs[i,i]   = 3/8*delta*homogeneous_fuel.absorption()
        if i==mesh_p-1:
            op_abs[i,i]   = 3/8*delta*homogeneous_fuel.absorption()
            op_abs[i,i-1] = 1/8*delta*homogeneous_fuel.absorption()

# Calculate the loss operator (absorption + leakage)

op_loss = op_abs + op_diff

# Initialize flux, eigenvalue and source

flux_distribution_iter = np.ones(mesh_p)
flux_acc = [flux_distribution_iter]

k_iter = 1
k_acc = [k_iter]

source = op_gain@flux_distribution_iter
source_acc = [source]
scaled_source = source/k_iter

# Solve the eigenvalue problem using the source iteration method

iteration = 0
max_iterations = 100
error = 0.0001
converged = False

while(not converged and max_iterations > iteration):
    iteration += 1
    
    flux_distribution_iter = np.linalg.solve(op_loss, scaled_source)
    flux_acc.append(flux_distribution_iter)
    
    source = op_gain@flux_distribution_iter
    
    k_iter = sum(source)/(sum(source_acc[-1])/k_acc[-1])
    source_acc.append(source)
    k_acc.append(k_iter)
    
    scaled_source = source/k_iter
    
    k_error = abs((k_acc[-1] - k_acc[-2])/k_acc[-1])
    
    if k_error < error: converged = True
    
# Plot results

sns.set()

plt.figure()
ax = plt.plot(k_acc)

plt.figure()
ax = plt.plot(flux_distribution_iter)

plt.show()
