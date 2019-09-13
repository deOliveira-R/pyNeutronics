#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 01:23:55 2019

@author: rodrigo
"""

from constants_n_factors import *
import isotopes as data
from macro import GroupMacro


class presetMaterial:
    def __init__(self, density, composition):
        self.density = density  # g/cm3
        self.composition = composition


UO2 = presetMaterial(10.4, {'U': 1, 'O': 2})
H2O = presetMaterial(0.75, {'H': 2, 'O': 1})


class Material:
    def __init__(self, density, fraction):
        self.composition = {}

        density_type = None
        if density is 'absolute':
            density_type = 'absolute' # fractions are absolute atom densities already
        elif density > 0:
            density_type = 'atom' # atom/cm3
        elif density < 0:
            density_type = 'mass' # g/cm3
        else:
            pass # TODO: implement an error catch

        # check all keys in fraction dictionary have same sign (atom or mass)
        test_key = list(fraction)[0] # fraction.keys()[0]
        atom_fraction = None
        if fraction[test_key] > 0:
            atom_fraction = True
        elif fraction[test_key] < 0:
            atom_fraction = False
        else:
            pass
            # TODO: implement error catch

        if atom_fraction:
            if any(fraction[key] < 0 for key in fraction):
                pass # TODO: implement error catch
            else:
                pass
        else:
            if any(fraction[key] > 0 for key in fraction):
                pass # TODO: implement error catch
            else:
                pass

        # if fractions are absolute, pass them to the composition dictionary
        if density_type is 'absolute':
            self.composition = fraction.copy()
        # if fractions are not absolute, calculate composition from density and fractions
        else:
            # After density and fraction classification, make everyone positive
            density = abs(density)
            fraction.update((key, abs(value)) for key, value in fraction.items())

            # Normalize fractions to 1
            factor = sum(fraction.values())
            fraction.update((key, value/factor) for key, value in fraction.items())

            # Calculate isotopic densities (isotope / cm3)
            if density_type is 'atom':
                if not atom_fraction:
                    self.massToAtom(fraction)

                self.composition = {key: value * density for key, value in fraction.items()}

            elif density_type is 'mass':
                if atom_fraction:
                    self.atomToMass(fraction)

                self.composition = {key: value * density * (Avogadro / getattr(data, key).mass)
                                    for key, value in fraction.items()}
            else:
                pass # TODO: wtf error


    @staticmethod
    def massToAtom(fraction):
        nf = lambda w, mm: w / mm
        total_den = sum(nf(value, getattr(data, key).mass) for key, value in fraction.items())
        fraction.update((key, nf(value, getattr(data, key).mass) / total_den) for key, value in fraction.items())

    @staticmethod
    def atomToMass(fraction):
        wf = lambda n, mm: n * mm
        total_den = sum(wf(value, getattr(data, key).mass) for key, value in fraction.items())
        fraction.update((key, wf(value, getattr(data, key).mass) / total_den) for key, value in fraction.items())

    @classmethod
    def mix(cls, mat1, mat2, ratio1to2):
        fraction2 = 1 / (ratio1to2 + 1)
        fraction1 = 1 - fraction2

        composition1 = {key: value * fraction1 for key, value in mat1.composition.items()}
        composition2 = {key: value * fraction2 for key, value in mat2.composition.items()}

        mixed = composition1.copy()
        for key in composition2:
            if key in mixed:
                mixed[key] += composition2[key]
            else:
                mixed[key] = composition2[key]

        return cls('absolute', mixed)


class NuclearMaterial:
    def __init__(self, lib, material):
        self._material = material
        self._partial_xs = {key: GroupMacro(lib, key, value) for key, value in self._material.composition.items()}

        self._scatter = sum(value.scatter() for value in self._partial_xs.values())
        self._capture = sum(value.capture() for value in self._partial_xs.values())
        self._fission = sum(value.fission() for value in self._partial_xs.values())
        self._absorption = sum(value.absorption() for value in self._partial_xs.values())
        self._total = sum(value.total() for value in self._partial_xs.values())
        self._nu_fission = sum(value.nu_fission() for value in self._partial_xs.values())

        self._diffusion = 1/(3*(self.total() - self.scatter() * self.scatter_angle()))

    def scatter_angle(self):
        angle = sum(value.scatter_angle() * value.scatter() for value in self._partial_xs.values())
        angle /= self.scatter()

        return angle

    def scatter(self):
        return self._scatter

    def capture(self):
        return self._capture

    def fission(self):
        return self._fission

    def absorption(self):
        return self._absorption

    def total(self):
        return self._total

    def nu_fission(self):
        return self._nu_fission

    def diffusion(self):
        return self._diffusion
