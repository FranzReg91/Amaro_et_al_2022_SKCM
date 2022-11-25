#!/usr/bin/env python
# coding: utf-8

# command example: python3 JSVD_Pfeffer_et_al_2019.py 4

import sys # add input varaibles from terminal (number of classes)
import math # to compute power of F. norm
import os # import operative system modules
import numpy as np
from pymanopt.manifolds import Product, Stiefel # compute product manifolds
from pymanopt import Problem
import autograd.scipy
from pymanopt.solvers import TrustRegions # load solver
import autograd.numpy as anp # to compute gradient and hessian automatically

os.chdir('../')

its=1000000
G = np.loadtxt("./results/dataSKCMjcmf_2D/preproc/RNA_prep.txt")
M = np.loadtxt("./results/dataSKCMjcmf_2D/preproc/Met_prep.txt")

I = np.shape(G)[0] # NB python is 0 indexed, matlab starts from one

# number of classes
R = int(sys.argv[1]) # number of classes
N1 = np.shape(G)[1]
N2 = np.shape(M)[1]

# math.pow is equal to 3^2 = 9
nG = math.pow(np.linalg.norm(G,'fro'), 2)
nM = math.pow(np.linalg.norm(M,'fro'), 2)

# create problem structure

manifold = Product([Stiefel(I,R), Stiefel(N1,R), Stiefel(N2,R)]) # U, V, W

# I'll use autograd to avoid Hessian computing
def cost(x):
    U=x[0]
    V=x[1]
    W=x[2]
    return anp.sum(-(anp.linalg.norm(anp.diag(anp.diag(anp.dot(anp.dot(anp.transpose(U),G),V))),'fro')**2)/nG -(anp.linalg.norm(anp.diag(anp.diag(anp.dot(anp.dot(anp.transpose(U),M),W))),'fro')**2)/nM)


problem = Problem(manifold=manifold, cost=cost, verbosity = 1)

# define solver

solver = TrustRegions(minstepsize=1e-16, mingradnorm=1e-12, maxiter=its, maxtime=float('inf')) 

# let Pymanopt do the rest
Xopt = solver.solve(problem)

# print U in text file
np.savetxt("./results/U_JSVD_Pfeffer_et_al_2019.txt", Xopt[0], delimiter = "\t")