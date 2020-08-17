#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 01:23:55 2019

@author: Rodrigo


"""

import isotopes

class GroupMicro:
    """
    Represents the concept of a microscopic XS of a SINGLE isotope
    for a discrete energy mesh.

    TODO: - Expand to include multiple energy groups
          - Add parametrization by background XS (for mixing) and temperature
    """

    def __init__(self, isotope, fyield, scatter, capture, fission, nubar):
        self.isotope = getattr(isotopes, isotope)
        self.mass = self.isotope.mass
        self.fyield = fyield

        self._scatter = scatter
        self._capture = capture
        self._fission = fission

        self._nubar = nubar

    def scatter(self):
        return self._scatter

    def capture(self):
        return self._capture

    def fission(self):
        return self._fission

    def nubar(self):
        return self._nubar

    #    def scatter(self, *args, energy=None, temperature, background):
    #        """
    #        Option of getting the XS for all energy values by giving only the
    #        background. If background is not give, it causes an error!
    #        """
    #        if energy is None:
    #            return self._scatter[background][temperature]
    #        else:
    #            return self._scatter[background][temperature][energy]

    #    def capture(self, *args, energy=None, temperature, background):
    #        """
    #        Option of getting the XS for all energy values by giving only the background
    #        If background is not give, it causes an error!
    #        """
    #        if energy is None:
    #            return self.capture[background][temperature]
    #        else:
    #            return self.capture[background][temperature][energy]

    #    def fission(self, *args, energy=None, temperature, background):
    #        """
    #        Option of getting the XS for all energy values by giving only the
    #        background. If background is not give, it causes an error!
    #        """
    #        if energy is None:
    #            return self._fission[background][temperature]
    #        else:
    #            return self._fission[background][temperature][energy]

    def absorption(self):
        return self.capture() + self.fission()

    def total(self):
        return self.absorption() + self.scatter()

    def nu_fission(self):
        return self.fission() * self.nubar()
