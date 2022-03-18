import matplotlib.pyplot as plt
import numpy as np

def plot_region(A=None, b=None, lb=None, ub=None, color='default', transparency=0.5, linetype='k', points=None, start_end_special=False):
    # plot_region plots a closed region or regions(s) in 2D/3D. the region(s) x is a subset of R2 or R3 s.t. Ax>=b
    # and lb <= x <= ub. A set of points can be included in the same plot if desired.

    #Inputs:
    # A: a matrix or set of matrices (if multiple, input as numpy arrays and concat along the last axis)
    # b, ub, lb: vectors or sets of vectors
    # color sets color ('r', 'g', etc); default indicates random
    # transparency between 0 and 1, default 0.5
    # points: set of points in the form of a matrix
    # start_end_special indicates whether to use a special marking for first/last point

    empty_vars = sum([A is None, b is None, lb is None, ub is None])
    if empty_vars > 2:
        raise Exception("Too few arguments for plot_region")

    # save dimensions to count number of affine eqn's
    if A is not None:
        if A.ndims > 2:
            m, n, p = A.shape()
        else:
            m, n = A.shape()
            p = 1

    # todo: dimensions check/fill missing b's etc

    # get number of bounds
    numb = 0
    if lb and lb.ndims > 1 and lb.shape()[-1] > numb:
        numb = lb.shape()[-1]
    if ub and ub.ndims > 1 and ub.shape()[-1] > numb:
        numb = ub.shape()[-1]

    # todo set colors for each eqn

    if n == 2:
        raise NotImplemented

    elif n==3:
        eq = [0 for i in range(p)]
        X = [0 for i in range(p)]

        # if p > 1:
        #     for pp in range(p):
        #         eq[pp] = np.zeros(3, 1)
        #         X[pp] = np.zeros(3, 1)
        #         Ap = A[:, :, pp]
        #         mp, np = Ap.shape
        #
        #         for i in range(mp-2):
        #             for j in range(i+1, mp-1):
        #                 for k in range(j+1, mp):


        




