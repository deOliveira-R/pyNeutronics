#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 02:01:07 2019

@author: Rodrigo

The purpose of this file is to generate a database of Isotopic data that can be
 shared by any microscopic cross-section library.
It is to be consulted when there is a need to inquire for data such as number
of nucleons, mass, decay constant, etc.

TODO: - Include more common data such as decay constant. Implement stuff from
       Nuclear Wallet Cards (https://www.nndc.bnl.gov/nudat2/indx_sigma.jsp)
      - Evaluate the idea of subclassing Isotope from Element and put common
       element data inside the common class somehow (natural isotopic
       composition, number of protons, name, etc).
      - Some helper function to calculate defect energy from conversion by
       fission or fusion.
      - Evaluate if moving the data here to a HDF file makes more sense.
"""

# class Element:
#     def __init__(self, name, symbol, protons):
#         self.name = name
#         self.symbol = symbol
#         self.protons = protons
#           isotopic composition


class Isotope:
    def __init__(self, name, symbol, protons, nucleons, mass):
        self.name = name
        self.symbol = symbol
        self.protons = protons
        self.nucleons = nucleons
        self.mass = mass

        # decay_constant = float
        # decay_chain = something...


H1    = Isotope('Hydrogen',  'H',  1,   1,   1.00782503223)
O16   = Isotope(  'Oxygen',  'O',  8,  16,  15.99491461957)
Xe135 = Isotope(   'Xenon', 'Xe', 54, 135, 134.9072278)
Sm149 = Isotope('Samarium', 'Sm', 62, 149, 148.9171921)
U235  = Isotope( 'Uranium',  'U', 92, 235, 235.0439301)
U238  = Isotope( 'Uranium',  'U', 92, 238, 238.0507884)
