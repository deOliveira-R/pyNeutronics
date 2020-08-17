#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 23:54:54 2019

@author: rodrigo
"""

from micro import GroupMicro

# Six Energy groups library for critical assemblies
# Isotopes: U233, U235, U238, Pu239 and Pu240

# Note that for 6 energy groups, 5 points between 0 and Infinity are needed to
# define 6 intervals
# Energy = [0.1, 0.4, 0.9, 1.4, 3.0] # in MeV


# Sixteen Energy groups library for critical assemblies
# Energy2 = [1E-7, 4E-7, 1E-6, 3E-6, 1E-5, 3E-5, 1E-4, 5.5E-4,
#            3E-3, 1.7E-2, 1E-1, 4E-1, 9E-1, 1.4, 3.0]

groups = 1

H1    = GroupMicro("H1",     # isotope
                   0,        # fission yield
                   3.3723,   # scatter
                   0,        # capture
                   0,        # fission
                   0)        # nubar

#                  isotope   yield   scatter      capture fission nubar
O16   = GroupMicro(  "O16",      0,  3.7156 ,      0     ,    0  ,    0)
Xe135 = GroupMicro("Xe135", 0.0640,  6.99340, 220000     ,    0  ,    0)
Sm149 = GroupMicro("Sm149", 0.0113,  8.65770,   7200     ,    0  ,    0)
U235  = GroupMicro( "U235",      0, 11.5860 ,     10.350 ,   47.5,    2.5)
U238  = GroupMicro( "U238",      0, 11.2221 ,      0.8900,    0  ,    0)
