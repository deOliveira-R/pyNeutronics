#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  19 13:58:00 2020

@author: Rodrigo

Implements 1 dimensional domain class to be used by solvers for geometrical description of the problem.
"""

from itertools import chain
import numpy as np
import unittest


class Domain:

    shapes = {'slab': ['slab', 'slb'],
              'cylindrical': ['cylindrical', 'cyl'],
              'spherical': ['spherical', 'sph']}

    cartesian = shapes['slab']
    curvilinear = shapes['cylindrical'] + shapes['spherical']

    def __init__(self, points, geometry: str):
        """
        Using an array of points and given shape, creates a domain by calculating:
         - delta_points: distance between points
         - centers: cell centers
         - delta_centers: distance between cell centers (useful to calculate finite difference derivatives at the interfaces)
         - surfaces: surface area at points between cells and at geometry boundaries
         - volumes: volume of cells

        :param points: array of points representing cell boundaries
        :param geometry: shape of domain
        """
        self.pre_checks(points, geometry)

        self.points = points
        self.delta_points = self.points[1:] - self.points[:-1]
        self.centers = self.points[:-1] + 0.5 * self.delta_points
        self.delta_centers = self.centers[1:] - self.centers[:-1]

        if geometry in self.shapes['slab']:
            # in slab it’s 1 everywhere
            self.surfaces = np.ones_like(points)
            # in slab its ri+1 - ri (which is delta at ri)
            self.volumes = self.delta_points
        elif geometry in self.shapes['cylindrical']:
            # in cylinder it is 2 pi r
            self.surfaces = 2.0 * np.pi * points
            # in cylinder its pi (r1^2 - r0^2)
            self.volumes = np.pi * (points[1:] ** 2 - points[:-1] ** 2)
        elif geometry in self.shapes['spherical']:
            # in sphere it is 4 pi r^2
            self.surfaces = 4.0 * np.pi * points ** 2
            # in sphere its 4/3 pi (ri+1^3 - ri^3)
            self.volumes = 4.0 / 3.0 * np.pi * (points[1:] ** 3 - points[:-1] ** 3)
        else:
            raise ValueError(f'Unspecified geometry type. Must be in {self.shapes}')

    @ classmethod
    def uniform(cls, length: float, intervals: int, geometry: str):
        """
        Helper initializer that simplifies creating a uniform domain of size "length" (starting from 0),
        composed of "intervals" cells and I+1 points (cell interfaces).

        :param length: size of domain (length of slab or radius of cylinder or sphere)
        :param intervals: number of cells to divide the domain
        :param geometry: shape of domain (i.e.: slab, cylinder or sphere)
        :return: class instance with uniform domain
        """

        delta = length / intervals
        points = np.arange(intervals + 1) * delta

        return cls(points, geometry)

    def pre_checks(self, points, geometry):
        assert self.strictly_increasing(points), \
            "Array of points must be strictly increasing."

        assert geometry in chain(*self.shapes.values()), \
            f"Geometry type '{geometry}' is unavailable. Available geometries are: {list(chain(*self.shapes.values()))}"

        if geometry in self.curvilinear:
            assert points[0] == 0.0,\
                f"First point of curvilinear geometries must be 0.0\n First point given is {points[0]}"

        return

    @staticmethod
    def strictly_increasing(arr) -> bool:
        return all(x < y for x, y in zip(arr, arr[1:]))


class DomainTest(unittest.TestCase):

    def setUp(self):
        self.length = 10
        self.intervals = 5
        self.uniform_slab = Domain.uniform(self.length, self.intervals, 'slb')
        self.uniform_cylinder = Domain.uniform(self.length, self.intervals, 'cyl')
        self.uniform_sphere = Domain.uniform(self.length, self.intervals, 'sph')

        self.nonuniform_grid = np.array([0.0, 1.0, 4.0, 5.5, 8.0, 10.0])
        self.nonuniform_intervals = (self.nonuniform_grid.size - 1)
        self.nonuniform_slab = Domain(self.nonuniform_grid, 'slb')
        self.nonuniform_cylinder = Domain(self.nonuniform_grid, 'cyl')
        self.nonuniform_sphere = Domain(self.nonuniform_grid, 'sph')

    def test_uniform_points(self):
        points = self.uniform_slab.points
        self.assertTrue(points.size == self.intervals + 1)
        reference_points = np.array([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
        self.assertTrue(np.allclose(points, reference_points))

    def test_uniform_delta_points(self):
        deltas = self.uniform_slab.delta_points
        self.assertTrue(deltas.size == self.intervals)
        reference_deltas = np.array([2.0, 2.0, 2.0, 2.0, 2.0])
        self.assertTrue(np.allclose(deltas, reference_deltas))

    def test_uniform_centers(self):
        centers = self.uniform_slab.centers
        self.assertTrue(centers.size == self.intervals)
        reference_points = np.array([1.0, 3.0, 5.0, 7.0, 9.0])
        self.assertTrue(np.allclose(centers, reference_points))

    def test_uniform_delta_centers(self):
        deltas = self.uniform_slab.delta_centers
        self.assertTrue(deltas.size == self.intervals - 1)
        reference_deltas = np.array([2.0, 2.0, 2.0, 2.0])
        self.assertTrue(np.allclose(deltas, reference_deltas))

    def test_uniform_surfaces_slab(self):
        slab_s = self.uniform_slab.surfaces
        reference_slab_s = np.ones_like(self.intervals + 1)
        self.assertTrue(np.allclose(slab_s, reference_slab_s))

    def test_uniform_volumes_slab(self):
        slab_v = self.uniform_slab.volumes
        reference_slab_v = np.full_like(self.intervals, self.length/self.intervals)
        self.assertTrue(np.allclose(slab_v, reference_slab_v))

    def test_uniform_surfaces_cylinder(self):
        cylinder_s = self.uniform_cylinder.surfaces
        reference_cylinder_s = np.array([0.0, 12.5664, 25.1327, 37.6991, 50.2655, 62.8319])
        self.assertTrue(np.allclose(cylinder_s, reference_cylinder_s))

    def test_uniform_volumes_cylinder(self):
        cylinder_v = self.uniform_cylinder.volumes
        reference_cylinder_v = np.array([12.5664, 37.6991, 62.8319, 87.9646, 113.097])
        self.assertTrue(np.allclose(cylinder_v, reference_cylinder_v))

    def test_uniform_surfaces_sphere(self):
        sphere_s = self.uniform_sphere.surfaces
        reference_sphere_s = np.array([0.0, 50.2655, 201.062, 452.389, 804.248, 1256.64])
        self.assertTrue(np.allclose(sphere_s, reference_sphere_s))

    def test_uniform_volumes_sphere(self):
        sphere_v = self.uniform_sphere.volumes
        reference_sphere_v = np.array([33.5103, 234.572, 636.696, 1239.88, 2044.13])
        self.assertTrue(np.allclose(sphere_v, reference_sphere_v))

    def test_nonuniform_delta_points(self):
        deltas = self.nonuniform_slab.delta_points
        self.assertTrue(deltas.size == self.nonuniform_intervals)
        reference_deltas = np.array([1.0, 3.0, 1.5, 2.5, 2.0])
        self.assertTrue(np.allclose(deltas, reference_deltas))

    def test_nonuniform_centers(self):
        centers = self.nonuniform_slab.centers
        self.assertTrue(centers.size == self.nonuniform_intervals)
        reference_centers = np.array([0.5, 2.5, 4.75, 6.75, 9.0])
        self.assertTrue(np.allclose(centers, reference_centers))

    def test_nonuniform_delta_centers(self):
        deltas = self.nonuniform_slab.delta_centers
        self.assertTrue(deltas.size == self.nonuniform_intervals - 1)
        reference_deltas = np.array([2.0, 2.25, 2.0, 2.25])
        self.assertTrue(np.allclose(deltas, reference_deltas))

    def test_nonuniform_surfaces_cylinder(self):
        cylinder_s = self.nonuniform_cylinder.surfaces
        reference_cylinder_s = np.array([0.0, 6.28319, 25.1327, 34.5575, 50.2655, 62.8319])
        self.assertTrue(np.allclose(cylinder_s, reference_cylinder_s))

    def test_nonuniform_volumes_cylinder(self):
        cylinder_v = self.nonuniform_cylinder.volumes
        reference_cylinder_v = np.array([3.14159, 47.1239, 44.7677, 106.029, 113.097])
        self.assertTrue(np.allclose(cylinder_v, reference_cylinder_v))

    def test_nonuniform_surfaces_sphere(self):
        sphere_s = self.nonuniform_sphere.surfaces
        reference_sphere_s = np.array([0.0, 12.5664, 201.062, 380.133, 804.248, 1256.637])
        self.assertTrue(np.allclose(sphere_s, reference_sphere_s))

    def test_nonuniform_volumes_sphere(self):
        sphere_v = self.nonuniform_sphere.volumes
        reference_sphere_v = np.array([4.18879, 263.894, 428.827, 1447.75, 2044.13])
        self.assertTrue(np.allclose(sphere_v, reference_sphere_v))


if __name__ == '__main__':
    unittest.main()
