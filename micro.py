import isotopes as data

class GroupMicro:
    """
    Represents the concept of a microscopic XS

    TODO: Parametrization by: energy group and background XS (for mixing)
    """

    def __init__(self, isotope, fyield, scatter, capture, fission, nubar):
        self.isotope = getattr(data, isotope)
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
