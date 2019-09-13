#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 02:01:07 2019

@author: rodrigo
"""

from constants_n_factors import *
import isotopes as data


# class Element:
#     def __init__(self, name, symbol, protons):
#         self.name = name
#         self.symbol = symbol
#         self.protons = protons
#
#
# class Isotope(Element):
#     def __init__(self, name, symbol, protons, neutrons, mass):
#         Element.__init__(name, symbol, protons)
#         self.neutrons = neutrons
#         self.mass = mass
#
#         # decay_constant = float
#         # decay_chain = something...
        

# class Material():
#     """
#     This should define a basic material, without nuclear data yet
#
#     It must have:
#         The composition of the material in stoichiometry in a dictionary
#
#     It must return:
#         The atom fraction for an isotope when given a string with its name
#         The mass fraction for an isotope when given a string with its name
#
#         Should support U as an argument or U235. For U, give total for all U isotopes
#
#         def split_symbol(s):
#             nucleon_number = s.lstrip(ascii_letters)
#             element = s[:len(s) - len(nucleon_number)]
#             return element, nucleon_number
#     """
#     def __init__(self, composition, mm, density):
#         self.composition = composition
#         self.mm = mm  # g/mol
#         self.density = density  # g/cm3
#
#     @classmethod
#     def mass_mixing(cls, mat1, mat2, ratio1to2):
#         mat2_fraction = 1/(ratio1to2 + 1)
#         mat1_fraction = 1 - mat2_fraction
#
#         density = mat1_fraction*mat1['density'] + mat2_fraction*mat2['density']
#
#         composition = {}
#         for i in mat1['composition']:
#             if i in composition:
#                 composition[i] += mat1_fraction*mat1['composition'][i]
#             else:
#                 composition[i] = mat1_fraction*mat1['composition'][i]
#
#         for i in mat2['composition']:
#             if i in composition:
#                 composition[i] += mat2_fraction*mat2['composition'][i]
#             else:
#                 composition[i] = mat2_fraction*mat2['composition'][i]
#
#         return cls(composition, mm, density)
#
#     def mol_density(self):
#         return self.density/self.mm  # mol/cm3
#
#     def molecular_density(self):
#         return self.mol_density()*avogadro  # molecules/cm3
    


    



    
#class nuclear_material():
#    def __init__(self, material):
        
        
        
        