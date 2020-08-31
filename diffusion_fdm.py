#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  19 13:58:00 2020

@author: Rodrigo

1 energy, 1 dimension problem. (mono-energetic slab)
Compute flux using a neutron diffusion approximation.

TODO: - Implement unit test with analytical solution
      - Implement time derivative using uniform mesh (linspace)
      - Test non-uniform time using diffusion stability condition
"""

import numpy as np
# import scipy.linalg as spl
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
# from scipy.special import j0
import unittest

from domain import Domain
import eigenvalue_solvers as es


class DomainOld:
    def __init__(self, length: float, intervals: int, geometry: str):
        """
        Create the points (cell edges) and centers for a uniform domain of size R and I cells

        :param length: size of domain
        :param intervals: number of cells
        """

        self.delta = float(length) / intervals
        self.centers = np.arange(intervals) * self.delta + 0.5 * self.delta
        self.points = np.arange(intervals + 1) * self.delta

        if geometry is 'slab' or geometry is 'slb':
            # in slab it’s 1 everywhere
            self.surfaces = np.ones_like(self.points)
            # #in slab its dr
            self.volumes = np.full_like(self.centers, self.delta)
        elif geometry is 'cylindrical' or geometry is 'cyl':
            # in cylinder it is 2 pi r
            self.surfaces = 2.0 * np.pi * self.points
            # in cylinder its pi (r^2 - r^2)
            self.volumes = np.pi * (self.points[1:] ** 2 - self.points[:-1] ** 2)
        elif geometry is 'spherical' or geometry is 'sph':
            # in sphere it is 4 pi r^2
            self.surfaces = 4.0 * np.pi * self.points ** 2
            # in sphere its 4/3 pi (r^3 - r^3)
            self.volumes = 4.0 / 3.0 * np.pi * (self.points[1:] ** 3 - self.points[:-1] ** 3)
        else:
            raise ValueError(
                'Unspecified geometry type. Must be:\n - slab (slb)\n - cylindrical (cyl)\n - spherical (sph)')


def mean_old(diffusion: np.ndarray, kind: str = 'harmonic') -> np.ndarray:
    """
    Finds the harmonic or arithmetic mean of diffusion constants at interfaces from cell values.

    :param diffusion: diffusion values in cells (size L)
    :param kind: selects between harmonic or arithmetic mean
    :return: array of mean values (size L - 1)
    """
    if kind is 'harmonic' or kind is 'har':
        diffusion_ = 2*diffusion[1:]*diffusion[:-1] / (diffusion[1:] + diffusion[:-1])
    elif kind is 'arithmetic' or kind is 'ari':
        diffusion_ = (diffusion[1:] + diffusion[:-1]) / 2
    else:
        raise ValueError(
            'Unspecified mean type. Must be:\n - harmonic (har)\n - arithmetic (ari)')

    return diffusion_


def DiffusionEigenvalueOriginalSymmetric(R: float, I: int, dif_func, sig_a_func, nu_fission_func, BC, geometry: str, epsilon: float = 1.0e-8):
    """
    Solve a neutron diffusion eigenvalue problem in a 1-D geometry using cell-averaged unknowns

    :param R: size of domain
    :param I: number of cells
    :param dif_func: name of function that returns diffusion coefficient for a given r
    :param sig_a_func: name of function that returns Sigma_a for a given r
    :param nu_fission_func: name of function that returns nu Sigma_f for a given r
    :param BC: Boundary Condition at r=R in form [A,B]
    :param geometry: shape of problem
    :param epsilon: tolerance
    :return:
            k: the multiplication factor of the system
            phi: the fundamental mode with norm 1
            centers: position at cell centers
    """
    # create the grid
    domain = DomainOld(R, I, geometry)
    delta = domain.delta
    centers = domain.centers
    S = domain.surfaces
    V = domain.volumes

    # This enforces Reflective BC without need to rely on accurate
    # evaluation of flux gradient at the inner boundary
    # Therefore it decreases the need for a finer mesh there
    S[0] = 0.0

    sig_a = sig_a_func(centers)
    nu_fission = nu_fission_func(centers)
    dif = dif_func(centers)

    D = mean_old(dif)
    D = extrapolate_to_boundaries(D)

    A = np.zeros((I + 1, I + 1))
    B = np.zeros_like(A)

    # Set up BC at R
    A[I, I] = (BC[0] * 0.5 + BC[1] / delta)
    A[I, I - 1] = (BC[0] * 0.5 - BC[1] / delta)

    # fill the matrix by running through the tri-diagonal
    for i, r in enumerate(centers):
        if i > 0:
            A[i, i - 1] = -S[i] * D[i] / (delta * V[i])
        A[i, i] = S[i] * D[i] / (delta * V[i]) \
                + S[i + 1] * D[i + 1] / (delta * V[i]) \
                + sig_a[i]
        A[i, i + 1] = -S[i + 1] * D[i + 1] / (delta * V[i])

        B[i, i] = nu_fission[i]

    # find eigenvalue
    l, phi = es.inverse_power(A, B, epsilon)

    k = 1.0 / l
    # remove last element of phi because it is outside the domain (ghost point)
    phi = phi[:I]
    return A, B, k, phi, centers


def extrapolate_to_boundaries(arr):
    arr = np.insert(arr, 0, arr[0])  # duplicate first element
    arr = np.append(arr, arr[-1])    # duplicate last element
    return arr


def mean(diffusion: np.ndarray, deltas: np.ndarray, kind: str = 'harmonic') -> np.ndarray:
    """
    Finds the harmonic or arithmetic mean of diffusion constants at interfaces from cell values using
    distance from cell centers to interfaces as weights.

    :param diffusion: diffusion values in cells (size L)
    :param deltas: cell lengths for spatially-weighted interpolation (size L)
    :param kind: selects between harmonic or arithmetic mean
    :return: array of mean values (size L - 1)

    TODO: need to check implementation of harmonic interpolation and implement unit test.
    """
    if kind is 'harmonic' or kind is 'har':
        diffusion_ = 1 / ((deltas[1:] / diffusion[1:] + deltas[:-1] / diffusion[:-1]) / (deltas[1:] + deltas[:-1]))
    elif kind is 'arithmetic' or kind is 'ari':
        diffusion_ = (diffusion[1:] * deltas[1:] + diffusion[:-1] * deltas[:-1]) / (deltas[1:] + deltas[:-1])
    else:
        raise ValueError(
            'Unspecified mean type. Must be:\n - harmonic (har)\n - arithmetic (ari)')

    return diffusion_


class DiffusionSource:
    """
    Solve a source driven neutron diffusion problem in a 1-D geometry using cell-averaged unknowns.
    """

    def __init__(self, domain: Domain, dif, sig_a, nu_fission, source, mean_kind: str = 'harmonic'):
        """
        Initialize solver from an instance of Domain, containing geometrical description of the problem
        and lists of nuclear data.

        :param domain: instance of Domain containing geometrical description
        :param dif: list of diffusion coefficients at cells
        :param sig_a: list of absorption cross-sections at cells
        :param nu_fission: list of neutron production cross-sections at cells
        :param source: list of neutron sources at cells
        :param mean_kind: type of mean to use for diffusion interpolation at cell boundaries (harmonic or arithmetic)
        """
        self.domain = domain

        # get cell values
        self.delta_p = self.domain.delta_points  # useful for interpolation of cell-defined values
        self.centers = self.domain.centers
        self.V = self.domain.volumes

        self.dif = dif
        self.sig_a = sig_a
        self.nu_fission = nu_fission
        self.source = source

        # get values at interface
        self.delta_c = self.domain.delta_centers  # distance between cell centers (useful for finite differences)
        self.S = self.domain.surfaces
        self.D = mean(dif, self.delta_p, kind=mean_kind)

        # These values are not defined at the outermost surface of the domain (boundaries)
        # By extrapolating delta_c, we created a "ghost cell" outside of the domain at the same distance from the
        # boundary as the boundary to the first cell center. The ghost cell approach is typical in finite differences
        # to apply boundary conditions when the cell values are defined at cell centers, not edges.
        # Extrapolation of diffusion coefficient is the current approach. Maybe it could be replaced by an effective
        # coefficient.
        self.delta_c = extrapolate_to_boundaries(self.delta_c)
        self.D = extrapolate_to_boundaries(self.D)

    @classmethod
    def from_position_function(cls, domain: Domain, dif_func, sig_a_func, nu_fission_func, source_func, mean_kind: str = 'harmonic'):
        """
        Helper function to initialize the solver from cross-sections given as functions that take an array of positions
        and return an array of cross-sections.

        :param domain: instance of Domain containing geometrical description
        :param dif_func: function that returns diffusion coefficients given positions
        :param sig_a_func: function that returns absorption cross-section given positions
        :param nu_fission_func: function that returns neutron production cross-section given positions
        :param source_func: function that returns source given positions
        :param mean_kind: type of mean to use for diffusion interpolation at cell boundaries (harmonic or arithmetic)
        :return: instance of DiffusionSource
        """
        dif = dif_func(domain.centers)
        sig_a = sig_a_func(domain.centers)
        nu_fission = nu_fission_func(domain.centers)
        source = source_func(domain.centers)

        return cls(domain, dif, sig_a, nu_fission, source, mean_kind=mean_kind)

    @classmethod
    def from_materials(cls, domain: Domain, materials, source_func, mean_kind: str = 'harmonic'):
        """
        Helper function to initialize the solver from list of Material instances to provide cross-sections.

        :param domain: instance of Domain containing geometrical description
        :param materials: list of instances of Material at each cell position
        :param source_func: function that returns source given positions
        :param mean_kind: type of mean to use for diffusion interpolation at cell boundaries (harmonic or arithmetic)
        :return: instance of DiffusionSource
        """
        dif = np.array([mat.diffusion() for mat in materials])
        sig_a = np.array([mat.absorption() for mat in materials])
        nu_fission = np.array([mat.nu_fission() for mat in materials])
        source = source_func(domain.centers)

        return cls(domain, dif, sig_a, nu_fission, source, mean_kind=mean_kind)

    # def calculate_current(self, phi):
    #     self.D[i] * (phi[i + 1] - phi[i]) / self.delta_c[i]

    def assemble(self, BC, BCO=(0, 1, 0)):
        """
        Assemble left hand side and right hand side matrices of the problem given its boundary conditions.

        :param BC: Boundary Condition at r=R in form (A, B, C)
        :param BCO: Boundary Condition at r=0 (origin). Defaults to reflective one.
                    Use (1, 0, 0) for flux 0 at boundary or
                    Use (0.24, D[0]/2, 0) for vacuum boundary (0 incoming current)
                    Use third value to impose source at boundaries
        :return:
                A: left hand side matrix containing neutron sinks
                b: right hand side vector containing neutron sources
        """

        # # create matrices with size of cell centers + 2 extras for boundary conditions
        # sz = len(self.centers)
        # A = np.zeros((sz + 2, sz + 2))
        # B = np.zeros_like(A)
        # b = np.zeros(sz + 2)
        #
        # # Set up boundary conditions at 0 and R on the first and last rows of matrix A
        # A[0, 0] = (BCO[0] * 0.5 + BCO[1] / self.delta_c[0])
        # A[0, 1] = (BCO[0] * 0.5 - BCO[1] / self.delta_c[0])
        # b[0] = BCO[2]
        # A[sz + 1, sz] = (0.5 * BC[0] - BC[1] / self.delta_c[-1])
        # A[sz + 1, sz + 1] = (0.5 * BC[0] + BC[1] / self.delta_c[-1])
        # b[sz+1] = BC[2]
        #
        # # fill the matrix by running through the tri-diagonal
        # # we avoid the first and last rows, since they are boundary conditions already set
        # for i, r in enumerate(self.centers):
        #     m = i + 1
        #
        #     backward = self.S[i] * self.D[i] / (self.delta_c[i] * self.V[i])
        #     forward = self.S[i + 1] * self.D[i + 1] / (self.delta_c[i + 1] * self.V[i])
        #
        #     A[m, m - 1] = -backward
        #     A[m, m] = backward + forward + self.sig_a[i] - self.nu_fission[i]
        #     A[m, m + 1] = -forward
        #
        #     B[m, m] = self.nu_fission[i]
        #     b[m] = self.source[i]

        # Set up boundary conditions at 0 and R on the first and last rows of matrix A
        bcoe = (BCO[0] * 0.5 + BCO[1] / self.delta_c[0])
        bco0 = (BCO[0] * 0.5 - BCO[1] / self.delta_c[0])
        bc0 = (0.5 * BC[0] - BC[1] / self.delta_c[-1])
        bce = (0.5 * BC[0] + BC[1] / self.delta_c[-1])

        backwards = self.S[:-1] * self.D[:-1] / (self.delta_c[:-1] * self.V)
        forwards = self.S[1:] * self.D[1:] / (self.delta_c[1:] * self.V)
        centrals = backwards + forwards + self.sig_a - self.nu_fission

        A = sps.diags([np.concatenate([-backwards, [bc0]]),
                        np.concatenate([[bcoe], centrals, [bce]]),
                        np.concatenate([[bco0], -forwards])],
                       [-1, 0, 1])

        b = np.concatenate([[BCO[2]], self.source, [BC[2]]])

        return A, b

    def solve(self, BC, BCO=(0, 1, 0), epsilon: float = 1.0e-8):
        """
        Solve source driven problem with given boundary condition and tolerance

        :param BC: Boundary Condition at r=R in form (A, B, C)
        :param BCO: Boundary Condition at r=0 (origin). Defaults to reflective one.
                    Use (1, 0, 0) for flux 0 at boundary or
                    Use (0.24, D[0]/2, 0) for vacuum boundary (0 incoming current)
                    Use third value to impose source at boundaries
        :param epsilon: tolerance on solution
        :return:
                phi: the fundamental mode with norm 1
                centers: position at cell centers
        TODO: implement selectable solver
        """
        A, b = self.assemble(BC, BCO)

        # solve problem
        # phi = np.linalg.solve(A, b)
        phi = spsl.spsolve(A.tocsr(), b)

        # remove first and last elements of phi because they are outside the domain (at ghost points)
        phi = phi[1:-1]
        return phi, self.centers


class DiffusionEigenvalue1E:
    """
    Solve a neutron diffusion eigenvalue problem in a 1-D geometry using cell-averaged unknowns.
    """

    def __init__(self, domain: Domain, dif, sig_a, nu_fission, mean_kind: str = 'harmonic'):
        """
        Initialize solver from an instance of Domain, containing geometrical description of the problem
        and lists of nuclear data.

        :param domain: instance of Domain containing geometrical description
        :param dif: list of diffusion coefficients at cells
        :param sig_a: list of absorption cross-sections at cells
        :param nu_fission: list of neutron production cross-sections at cells
        :param mean_kind: type of mean to use for diffusion interpolation at cell boundaries (harmonic or arithmetic)
        """
        self.domain = domain

        # get cell values
        self.delta_p = self.domain.delta_points  # useful for interpolation of cell-defined values
        self.centers = self.domain.centers
        self.V = self.domain.volumes

        self.dif = dif
        self.sig_a = sig_a
        self.nu_fission = nu_fission

        # get values at interface
        self.delta_c = self.domain.delta_centers  # distance between cell centers (useful for finite differences)
        self.S = self.domain.surfaces
        self.D = mean(dif, self.delta_p, kind=mean_kind)

        # These values are not defined at the outermost surface of the domain (boundaries)
        # By extrapolating delta_c, we created a "ghost cell" outside of the domain at the same distance from the
        # boundary as the boundary to the first cell center. The ghost cell approach is typical in finite differences
        # to apply boundary conditions when the cell values are defined at cell centers, not edges.
        # Extrapolation of diffusion coefficient is the current approach. Maybe it could be replaced by an effective
        # coefficient.
        self.delta_c = extrapolate_to_boundaries(self.delta_c)
        self.D = extrapolate_to_boundaries(self.D)

    @classmethod
    def from_position_function(cls, domain: Domain, dif_func, sig_a_func, nu_fission_func, mean_kind: str = 'harmonic'):
        """
        Helper function to initialize the solver from cross-sections given as functions that take an array of positions
        and return an array of cross-sections.

        :param domain: instance of Domain containing geometrical description
        :param dif_func: function that returns diffusion coefficients given positions
        :param sig_a_func: function that returns absorption cross-section given positions
        :param nu_fission_func: function that returns neutron production cross-section given positions
        :param mean_kind: type of mean to use for diffusion interpolation at cell boundaries (harmonic or arithmetic)
        :return: instance of DiffusionEigenvalue
        """
        sig_a = sig_a_func(domain.centers)
        nu_fission = nu_fission_func(domain.centers)
        dif = dif_func(domain.centers)

        return cls(domain, dif, sig_a, nu_fission, mean_kind=mean_kind)

    @classmethod
    def from_materials(cls, domain: Domain, materials, mean_kind: str = 'harmonic'):
        """
        Helper function to initialize the solver from list of Material instances to provide cross-sections.

        :param domain: instance of Domain containing geometrical description
        :param materials: list of instances of Material at each cell position
        :param mean_kind: type of mean to use for diffusion interpolation at cell boundaries (harmonic or arithmetic)
        :return: instance of DiffusionEigenvalue
        """
        sig_a = np.array([mat.absorption() for mat in materials])
        nu_fission = np.array([mat.nu_fission() for mat in materials])
        dif = np.array([mat.diffusion() for mat in materials])

        return cls(domain, dif, sig_a, nu_fission, mean_kind=mean_kind)

    # def calculate_current(self, phi):
    #     self.D[i] * (phi[i + 1] - phi[i]) / self.delta_c[i]

    def assemble(self, BC, BCO=(0, 1)):
        """
        Assemble left hand side and right hand side matrices of the problem given its boundary conditions.

        :param BC: Boundary Condition at r=R in form [A,B]
        :param BCO: Boundary Condition at r=0 (origin). Defaults to reflective one.
                    Use (1, 0) for flux 0 at boundary or
                    Use (0.24, D[0]/2) for vacuum boundary (0 incoming current)
        :return:
                A: left hand side matrix containing neutron sinks
                B: right hand side matrix containing neutron sources
        """

        # # create matrices with size of cell centers + 2 extras for boundary conditions
        # sz = len(self.centers)
        # A = np.zeros((sz + 2, sz + 2))
        # B = np.zeros_like(A)
        #
        # # Set up boundary conditions at 0 and R on the first and last rows of matrix A
        #
        # A[0, 0] = (BCO[0] * 0.5 + BCO[1] / self.delta_c[0])
        # A[0, 1] = (BCO[0] * 0.5 - BCO[1] / self.delta_c[0])
        # A[sz + 1, sz] = (0.5 * BC[0] - BC[1] / self.delta_c[-1])
        # A[sz + 1, sz + 1] = (0.5 * BC[0] + BC[1] / self.delta_c[-1])
        #
        # # fill the matrix by running through the tri-diagonal
        # # we avoid the first and last rows, since they are boundary conditions already set
        # for i, r in enumerate(self.centers):
        #     m = i + 1
        #
        #     backward = self.S[i] * self.D[i] / (self.delta_c[i] * self.V[i])
        #     forward = self.S[i + 1] * self.D[i + 1] / (self.delta_c[i + 1] * self.V[i])
        #
        #     A[m, m - 1] = -backward
        #     A[m, m] = backward + forward + self.sig_a[i]
        #     A[m, m + 1] = -forward
        #
        #     B[m, m] = self.nu_fission[i]

        bcoe = (BCO[0] * 0.5 + BCO[1] / self.delta_c[0])
        bco0 = (BCO[0] * 0.5 - BCO[1] / self.delta_c[0])
        bc0 = (0.5 * BC[0] - BC[1] / self.delta_c[-1])
        bce = (0.5 * BC[0] + BC[1] / self.delta_c[-1])

        backwards = self.S[:-1] * self.D[:-1] / (self.delta_c[:-1] * self.V)
        forwards = self.S[1:] * self.D[1:] / (self.delta_c[1:] * self.V)
        centrals = backwards + forwards + self.sig_a

        A = sps.diags([np.concatenate([-backwards, [bc0]]),
                       np.concatenate([[bcoe], centrals, [bce]]),
                       np.concatenate([[bco0], -forwards])],
                      [-1, 0, 1])

        B = sps.diags(np.concatenate([[0.0], self.nu_fission, [0.0]]))

        return A, B

    def solve(self, BC, BCO=(0, 1), normalization_power=None, epsilon: float = 1.0e-8):
        """
        Solve eigenvalue problem with given boundary condition and tolerance

        :param BC: Boundary Condition at r=R in form [A,B]
        :param BCO: Boundary Condition at r=0 (origin). Defaults to reflective one.
                    Use (1, 0) for flux 0 at boundary or
                    Use (0.24, D[0]/2) for vacuum boundary (0 incoming current)
        :param epsilon: tolerance on eigenvalue
        :return:
                k: the multiplication factor of the system
                phi: the fundamental mode with norm 1
                centers: position at cell centers
        TODO: implement selectable solver (give an eigenvalue problem solver as argument?)
        """
        A, B = self.assemble(BC, BCO)

        # find eigenvalue
        # l, phi = spl.eig(A, B)

        # l, phi = es.inverse_power(A, B, epsilon)
        # l, phi = es.inverse_power_lu(A, B, epsilon)
        # l, phi = es.inverse_power_linear(A, B, epsilon)

        # l, phi = spl.eig(A.toarray(), B.toarray())
        # l, phi = es.inverse_power(A.tocsr(), B, epsilon)
        # l, phi = es.inverse_power_lu(A.tocsc(), B, epsilon)
        l, phi = es.inverse_power_linear(A, B, epsilon)

        # Need to add some normalization with physical sense here.
        # Either generate a distribution that sums to 1
        # phi /= phi.sum()
        # Or use phi, kappa fission and power to get the absolute flux

        k = 1.0 / l
        # remove first and last elements of phi because they are outside the domain (at ghost points)
        phi = phi[1:-1]
        return k, phi, self.centers


class DiffusionEigenvalueMG:
    """
    Solve a neutron diffusion eigenvalue problem in a 1-D geometry using cell-averaged unknowns.
    """

    def __init__(self, domain: Domain, dif, rem, nu_fission, scatmat, mean_kind: str = 'harmonic'):
        """
        Initialize solver from an instance of Domain, containing geometrical description of the problem
        and lists of nuclear data.

        :param domain: instance of Domain containing geometrical description
        :param dif: list of diffusion coefficients at cells
        :param sig_a: list of absorption cross-sections at cells
        :param nu_fission: list of neutron production cross-sections at cells
        :param mean_kind: type of mean to use for diffusion interpolation at cell boundaries (harmonic or arithmetic)
        """
        self.domain = domain

        # get cell values
        self.delta_p = self.domain.delta_points  # useful for interpolation of cell-defined values
        self.centers = self.domain.centers
        self.V = self.domain.volumes

        self.dif = dif
        self.rem = rem
        self.nu_fission = nu_fission
        self.scatmat = scatmat

        # get values at interface
        self.delta_c = self.domain.delta_centers  # distance between cell centers (useful for finite differences)
        self.S = self.domain.surfaces
        self.D = mean(dif, self.delta_p, kind=mean_kind)

        # These values are not defined at the outermost surface of the domain (boundaries)
        # By extrapolating delta_c, we created a "ghost cell" outside of the domain at the same distance from the
        # boundary as the boundary to the first cell center. The ghost cell approach is typical in finite differences
        # to apply boundary conditions when the cell values are defined at cell centers, not edges.
        # Extrapolation of diffusion coefficient is the current approach. Maybe it could be replaced by an effective
        # coefficient.
        self.delta_c = extrapolate_to_boundaries(self.delta_c)
        self.D = extrapolate_to_boundaries(self.D)

    @classmethod
    def from_position_function(cls, domain: Domain, dif_func, rem_func, nu_fission_func, scatmat_func, mean_kind: str = 'harmonic'):
        """
        Helper function to initialize the solver from cross-sections given as functions that take an array of positions
        and return an array of cross-sections.

        :param domain: instance of Domain containing geometrical description
        :param dif_func: function that returns diffusion coefficients given positions
        :param sig_a_func: function that returns absorption cross-section given positions
        :param nu_fission_func: function that returns neutron production cross-section given positions
        :param mean_kind: type of mean to use for diffusion interpolation at cell boundaries (harmonic or arithmetic)
        :return: instance of DiffusionEigenvalue
        """
        dif = dif_func(domain.centers)
        rem = rem_func(domain.centers)
        nu_fission = nu_fission_func(domain.centers)
        scatmat = scatmat_func(domain.centers)

        return cls(domain, dif, rem, nu_fission, scatmat, mean_kind=mean_kind)

    @classmethod
    def from_materials(cls, domain: Domain, materials, mean_kind: str = 'harmonic'):
        """
        Helper function to initialize the solver from list of Material instances to provide cross-sections.

        :param domain: instance of Domain containing geometrical description
        :param materials: list of instances of Material at each cell position
        :param mean_kind: type of mean to use for diffusion interpolation at cell boundaries (harmonic or arithmetic)
        :return: instance of DiffusionEigenvalue
        """
        dif = np.array([mat.diffusion() for mat in materials])
        rem = np.array([mat.absorption() for mat in materials])
        nu_fission = np.array([mat.nu_fission() for mat in materials])
        scatmat = np.array([mat.scattering_matrix() for mat in materials])


        return cls(domain, dif, rem, nu_fission, scatmat, mean_kind=mean_kind)

    # def calculate_current(self, phi):
    #     self.D[i] * (phi[i + 1] - phi[i]) / self.delta_c[i]

    def assemble(self, BC, BCO=(0, 1)):
        """
        Assemble left hand side and right hand side matrices of the problem given its boundary conditions.

        :param BC: Boundary Condition at r=R in form [A,B]
        :param BCO: Boundary Condition at r=0 (origin). Defaults to reflective one.
                    Use (1, 0) for flux 0 at boundary or
                    Use (0.24, D[0]/2) for vacuum boundary (0 incoming current)
        :return:
                A: left hand side matrix containing neutron sinks
                B: right hand side matrix containing neutron sources
        """

        # create matrices with size of cell centers + 2 extras for boundary conditions
        sz = len(self.centers)
        A = np.zeros((sz + 2, sz + 2))
        B = np.zeros_like(A)

        # Set up boundary conditions at 0 and R on the first and last rows of matrix A
        A[0, 0] = (BCO[0] * 0.5 + BCO[1] / self.delta_c[0])
        A[0, 1] = (BCO[0] * 0.5 - BCO[1] / self.delta_c[0])
        A[sz + 1, sz] = (0.5 * BC[0] - BC[1] / self.delta_c[-1])
        A[sz + 1, sz + 1] = (0.5 * BC[0] + BC[1] / self.delta_c[-1])

        # fill the matrix by running through the tri-diagonal
        # we avoid the first and last rows, since they are boundary conditions already set
        for i, r in enumerate(self.centers):
            m = i + 1

            backward = self.S[i] * self.D[i] / (self.delta_c[i] * self.V[i])
            forward = self.S[i + 1] * self.D[i + 1] / (self.delta_c[i + 1] * self.V[i])

            A[m, m - 1] = -backward
            A[m, m] = backward + forward + self.rem[i]
            A[m, m + 1] = -forward

            B[m, m] = self.nu_fission[i]

        return A, B

    def solve(self, BC, BCO=(0, 1), normalization_power=None, epsilon: float = 1.0e-8):
        """
        Solve eigenvalue problem with given boundary condition and tolerance

        :param BC: Boundary Condition at r=R in form [A,B]
        :param BCO: Boundary Condition at r=0 (origin). Defaults to reflective one.
                    Use (1, 0) for flux 0 at boundary or
                    Use (0.24, D[0]/2) for vacuum boundary (0 incoming current)
        :param epsilon: tolerance on eigenvalue
        :return:
                k: the multiplication factor of the system
                phi: the fundamental mode with norm 1
                centers: position at cell centers
        TODO: implement selectable solver (give an eigenvalue problem solver as argument?)
        """
        A, B = self.assemble(BC, BCO)

        # find eigenvalue
        # l, phi = sp.linalg.eig(A, B)
        l, phi = es.inverse_power(A, B, epsilon)
        # l, phi = inverse_power_lu(A, B, epsilon)
        # l, phi = inverse_power_bicgstab(A, B, epsilon)

        # Need to add some normalization with physical sense here.
        # Either generate a distribution that sums to 1
        # phi /= phi.sum()
        # Or use phi, kappa fission and power to get the absolute flux

        k = 1.0 / l
        # remove first and last elements of phi because they are outside the domain (at ghost points)
        phi = phi[1:-1]
        return k, phi, self.centers


class FunctionTest(unittest.TestCase):

    def setUp(self):
        self.length = 10
        self.intervals = 5
        self.uniform_slab = Domain.uniform(self.length, self.intervals, 'slb')

        self.nonuniform_grid = np.array([0.0, 1.0, 4.0, 5.5, 8.0, 10.0])
        self.nonuniform_intervals = (self.nonuniform_grid.size - 1)
        self.nonuniform_slab = Domain(self.nonuniform_grid, 'slb')

    def test_uniform_delta_centers(self):
        deltas = self.uniform_slab.delta_centers
        self.assertTrue(deltas.size == self.intervals - 1)
        reference_deltas = np.array([2.0, 2.0, 2.0, 2.0])
        self.assertTrue(np.allclose(deltas, reference_deltas))
        extrapolated_deltas = extrapolate_to_boundaries(deltas)
        self.assertTrue(extrapolated_deltas.size == self.intervals + 1)
        reference_deltas = np.array([2.0, 2.0, 2.0, 2.0, 2.0, 2.0])
        self.assertTrue(np.allclose(extrapolated_deltas, reference_deltas))

    def test_nonuniform_delta_centers(self):
        deltas = self.nonuniform_slab.delta_centers
        self.assertTrue(deltas.size == self.nonuniform_intervals - 1)
        reference_deltas = np.array([2.0, 2.25, 2.0, 2.25])
        self.assertTrue(np.allclose(deltas, reference_deltas))
        extrapolated_deltas = extrapolate_to_boundaries(deltas)
        self.assertTrue(extrapolated_deltas.size == self.nonuniform_intervals + 1)
        reference_deltas = np.array([2.0, 2.0, 2.25, 2.0, 2.25, 2.25])
        self.assertTrue(np.allclose(extrapolated_deltas, reference_deltas))


# class DiffusionTest(unittest.TestCase):
#
#     def setUp(self):
#         self.length = 10
#         self.intervals = 5
#         self.uniform_slab = Domain.uniform(self.length, self.intervals, 'slb')
#
#         self.nonuniform_grid = np.array([0.0, 1.0, 4.0, 5.5, 8.0, 10.0])
#         self.nonuniform_intervals = (self.nonuniform_grid.size - 1)
#         self.nonuniform_slab = Domain(self.nonuniform_grid, 'slb')
#
#     def test_uniform_points(self):
#         points = self.uniform_slab.points
#         self.assertTrue(points.size == self.intervals + 1)
#         reference_points = np.array([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
#         self.assertTrue(np.allclose(points, reference_points))


if __name__ == '__main__':
    unittest.main()
