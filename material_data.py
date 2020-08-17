#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 01:23:55 2019

@author: Rodrigo

TODO: - This REALLY needs unit testing! Specially the mixing.
"""

from constants_n_factors import *
import isotopes
from macro import GroupMacro

import numpy as np
import unittest


class PresetMaterial:
    def __init__(self, density, composition):
        self.density = density  # g/cm3
        self.composition = composition


UO2 = PresetMaterial(10.4, {'U': 1, 'O': 2})
H2O = PresetMaterial(0.75, {'H': 2, 'O': 1})


class Material:
    def __init__(self, density, fraction: dict):
        self.composition = {}

        # check is density is a string or a float of negative or positive value
        # to decide how to dealt with the fraction dictionary
        density_type = None
        if density is 'absolute':
            density_type = 'absolute'  # fractions are absolute atom densities already
        elif density > 0:
            density_type = 'atom'  # atom/cm3
        elif density < 0:
            density_type = 'mass'  # g/cm3
        else:
            pass  # TODO: implement an error catch

        # get first key as a reference to see if all other keys have the same sign
        test_key = list(fraction)[0]  # fraction.keys()[0]
        atom_fraction = None
        if fraction[test_key] > 0:
            atom_fraction = True
        elif fraction[test_key] < 0:
            atom_fraction = False
        else:
            pass
            # TODO: implement error catch

        # check all keys in fraction dictionary have same sign as the first key
        # because it's not possible to deal with a mix of atom and mass fractions
        if atom_fraction:
            if any(fraction[key] < 0 for key in fraction):
                pass  # TODO: implement error catch
            else:
                pass
        else:
            if any(fraction[key] > 0 for key in fraction):
                pass  # TODO: implement error catch
            else:
                pass

        # if fractions are absolute, pass them to the composition dictionary
        if density_type is 'absolute':
            self.composition = fraction.copy()
        # if fractions are not absolute, calculate absolute composition
        # from density and fractions
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
                    self.mass_to_atom(fraction)

                self.composition = {key: value * density for key, value in fraction.items()}

            elif density_type is 'mass':
                if atom_fraction:
                    self.atom_to_mass(fraction)

                self.composition = {key: value * density * (Avogadro / getattr(isotopes, key).mass)
                                    for key, value in fraction.items()}
            else:
                pass  # TODO: wtf error

    @staticmethod
    def convert(from_to: str, fraction: dict):
        if from_to is 'massToAtom':
            factor = lambda w, mm: w / mm
        elif from_to is 'atomToMass':
            factor = lambda n, mm: n * mm
        else:
            pass  # TODO: error related to unexpected fromTo value

        total_density = sum(factor(value, getattr(isotopes, key).mass) for key, value in fraction.items())
        fraction.update((key, factor(value, getattr(isotopes, key).mass) / total_density)
                        for key, value in fraction.items())

    @staticmethod
    def mass_to_atom(fraction: dict):
        nf = lambda w, mm: w / mm
        total_den = sum(nf(value, getattr(isotopes, key).mass) for key, value in fraction.items())
        fraction.update((key, nf(value, getattr(isotopes, key).mass) / total_den) for key, value in fraction.items())

    @staticmethod
    def atom_to_mass(fraction: dict):
        wf = lambda n, mm: n * mm
        total_den = sum(wf(value, getattr(isotopes, key).mass) for key, value in fraction.items())
        fraction.update((key, wf(value, getattr(isotopes, key).mass) / total_den) for key, value in fraction.items())

    @classmethod
    def mix(cls, material1, material2, ratio1to2: float):
        """
        This method receives 2 Material objects and a mix ratio and returns another Material object from the mix.

        :param material1: Material
        :param material2: Material
        :param ratio1to2: float
        :return: Material
        """

        # calculate fractions from ratios
        fraction2 = 1 / (ratio1to2 + 1)
        fraction1 = 1 - fraction2

        # scale compositions by the calculated fractions, since the composition dictionary
        # is already in isotope atomic density.
        composition1 = {key: value * fraction1 for key, value in material1.composition.items()}
        composition2 = {key: value * fraction2 for key, value in material2.composition.items()}

        # add the scaled compositions in the mixed dictionary
        mixed = composition1.copy()
        for key in composition2:
            if key in mixed:
                mixed[key] += composition2[key]
            else:
                mixed[key] = composition2[key]

        # return the new mixed material by calling the initializer of this Material class
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


class MaterialTest(unittest.TestCase):

    def setUp(self):
        self.fraction_atom = {'H1': 2, 'O16': 1}
        self.fraction_atom_normalized = {'H1': 0.6666666666666666, 'O16': 0.3333333333333333}
        self.fraction_mass_normalized = {'H1': 0.11191487328808075, 'O16': 0.8880851267119192}

    def test_atomToMass(self):
        to_mass = self.fraction_atom.copy()
        Material.atom_to_mass(to_mass)
        self.assertTrue(np.allclose(list(to_mass.values()),
                                    list(self.fraction_mass_normalized.values()),
                                    rtol=1e-5))

    def test_massToAtom(self):
        to_atom = self.fraction_mass_normalized.copy()
        Material.mass_to_atom(to_atom)
        self.assertTrue(np.allclose(list(to_atom.values()),
                                    list(self.fraction_atom_normalized.values()),
                                    rtol=1e-5))

    def test_conversionHypothesis(self):
        hypothesis = self.fraction_atom.copy()
        Material.atom_to_mass(hypothesis)
        Material.mass_to_atom(hypothesis)
        self.assertTrue(np.allclose(list(hypothesis.values()),
                                    list(self.fraction_atom_normalized.values()),
                                    rtol=1e-5))


if __name__ == '__main__':
    unittest.main()
