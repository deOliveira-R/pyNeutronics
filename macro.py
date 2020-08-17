#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 01:23:55 2019

@author: Rodrigo


"""

from constants_n_factors import barn_to_cm2

class GroupMacro:
    """
    Represents the concept of macroscopic XS, or the attenuation caused by a
    nucleon upon an impinging beam.
    """

    def __init__(self, lib, isotope, density):
        self.micro = getattr(lib, isotope)
        self.density = density

        self._scatter = self.micro.scatter() * density * barn_to_cm2
        self._capture = self.micro.capture() * density * barn_to_cm2
        self._fission = self.micro.fission() * density * barn_to_cm2
        self._nu_fission = self.micro.nu_fission() * density * barn_to_cm2

        self._scatter_angle = 2/(3*self.micro.mass)

    def update(self):
        pass

    def scatter(self):
        return self._scatter

    def capture(self):
        return self._capture

    def fission(self):
        return self._fission

    def absorption(self):
        return self.capture() + self.fission()

    def total(self):
        return self.absorption() + self.scatter()

    def nu_fission(self):
        return self._nu_fission

    def scatter_angle(self):
        return self._scatter_angle

    # def diffusion(self):
    #     return 1 / (3 * (self.total() - self.scatter * self.scatter_angle()))
