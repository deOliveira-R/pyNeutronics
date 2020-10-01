#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 19:51:40 2020

@author: Rodrigo
"""
import numpy as np
import scipy.linalg as spl
import scipy.sparse.linalg as sps
import unittest


def inverse_power_lu(A, B, epsilon=1.0e-6, debug=False):
    """
    Solve the generalized eigenvalue problem
    Ax = l B x using inverse power iteration Inputs

    :param A: The LHS matrix (must be invertible)
    :param B: The RHS matrix
    :param epsilon: tolerance on eigenvalue
    :param debug: print warnings
    :return:
            l: the smallest eigenvalue of the problem
            x: the associated eigenvector
    """
    N, M = A.shape
    assert N == M
    # generate initial guess
    x = np.random.random(N)
    # x = np.ones((N))
    x = x / spl.norm(x)  # make norm(x)==1

    # compute LU factorization of A
    # lu, piv = spl.lu_factor(A)
    solve = sps.factorized(A)

    x_old = x
    l_old = 0
    converged = False

    iteration = 1
    while not converged:
        # x = spl.lu_solve((lu, piv), np.dot(B, x_old))
        x = solve(B.dot(x_old))
        l = spl.norm(x)
        sign = x[0] / x_old[0] / l
        x_old = x / l

        converged = (np.fabs(l - l_old) < epsilon)
        l_old = l
        if debug:
            print("Iteration:", iteration, "\tMagnitude of l =", 1.0/l)
        iteration += 1

    return sign/l, x


def inverse_power_linear(A, B, epsilon=1.0e-6, solver=sps.bicgstab, debug=False):
    """
    Solve the generalized eigenvalue problem
    Ax = l B x using inverse power iteration Inputs

    :param A: The LHS matrix (must be invertible)
    :param B: The RHS matrix
    :param epsilon: tolerance on eigenvalue
    :param solver: sparse linear solver to use for the solution
    :param debug: print warnings
    :return:
            l: the smallest eigenvalue of the problem
            x: the associated eigenvector
    """
    N, M = A.shape
    assert N == M
    # generate initial guess
    x = np.random.random(N)
    # x = np.ones((N))
    x = x / spl.norm(x)  # make norm(x)==1

    x_old = x
    l_old = 1
    converged = False

    iteration = 1
    while not converged:
        x = solver(A, B.dot(x_old))

        l = spl.norm(x[0])
        sign = x[0][0] / x_old[0] / l
        x_old = x[0] / l

        converged = (np.fabs(l - l_old) < epsilon)
        l_old = l
        if debug:
            print("Iteration:", iteration, "\tMagnitude of l =", 1.0/l)
        iteration += 1

    return sign/l, x[0]


def inverse_power(A, B, epsilon=1.0e-6, debug=False):
    """
    Solve the generalized eigenvalue problem
    Ax = l B x using inverse power iteration Inputs

    :param A: The LHS matrix (must be invertible)
    :param B: The RHS matrix
    :param epsilon: tolerance on eigenvalue
    :param debug: print warnings
    :return:
            l: the smallest eigenvalue of the problem
            x: the associated eigenvector
    """
    N, M = A.shape
    assert N == M
    # generate initial guess
    x = np.random.random(N)
    # x = np.ones((N))
    x = x / spl.norm(x)  # make norm(x)==1

    x_old = x
    l_old = 1
    converged = False

    iteration = 1
    while not converged:
        # x = spl.solve(A, np.dot(B, x_old))
        x = sps.spsolve(A, B.dot(x_old))
        l = spl.norm(x)
        sign = x[0] / x_old[0] / l
        x_old = x / l

        converged = (np.fabs(l - l_old) < epsilon)
        l_old = l
        if debug:
            print("Iteration:", iteration, "\tMagnitude of l =", 1.0/l)
        iteration += 1

    return sign/l, x


# class SolverTest(unittest.TestCase):
#
#     def setUp(self):
#         self.lhs = np.identity(2)
#         self.lhs[1, 1] = 0.1
#         self.rhs = np.identity(2)
#         self.l_expected = 0.1
#
#     def test_inversePower(self):
#         l, x = inverse_power(self.lhs, self.rhs)
#         print(l, x)
#         self.assertAlmostEqual(l, self.l_expected, delta=1E-6)


if __name__ == '__main__':
    unittest.main()
