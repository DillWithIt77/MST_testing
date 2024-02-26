import numpy as np
import math
import sympy
from itertools import product, combinations
import itertools
import networkx as nx
import matplotlib.pyplot as plt

import sympy
import cdd
import time

#####chase's naive alrogithm circuit code
#naive approach for enumerating C(A,B)
#Input: A_eq and B_ineq are (m_a x n) and (m x n) numpy arrays
#Output: list of circuits in C(A,B) given by n-dimensional numpy arrays
def enumerate_circuits(B_ineq, A_eq=None):
    B = sympy.Matrix(B_ineq)
    m,n = B.shape
    # print(B)
    r = 0
    if A_eq is not None:
        A = sympy.Matrix(A_eq)
        A, pivot_columns = A.rref()    #use reduced echelon form of A
        r = len(pivot_columns)         #r is the rank of A

    circuits = []
    pos_circs = []
    for I in itertools.combinations(range(m),n-r-1):
        B_I = B[I,:]
        
        if A_eq is not None:
            D = A.col_join(B_I)
        else:
            D = B_I  
            
        ker_D = D.nullspace()
        if len(ker_D) == 1:   #circuit direction is found iff null space of D is one-dimensional
            # print('B_I rows',B_I)
            g = np.array(ker_D[0]).flatten()
            ###I think the .q implies the matrices need to come in integers/rationals? 
            ###It got upset when I had a float as an entry
            g = g*sympy.lcm([g[i].q for i in range(n) if g[i] != 0]) #normalize to coprime integers
            # print(f'circuit: {g}')
            g_is_duplicate = False
            for y in circuits:
                if np.array_equal(y, g):
                    g_is_duplicate = True
            if not g_is_duplicate:
                circuits.append(g)
                pos_circs.append(g)
                # print('circuit was added')
                ###commented out since don't necessarily care about having both directions
                circuits.append(-1*g)
                
    return (np.array(circuits),np.array(pos_circs))


def determine_y_vars(m,sign=None):
    y_vars = []
    if sign is None:
        for i in range(m):
            y_vars.append([i, 1])
            y_vars.append([i, -1])
    else:
        for i in range(m):
            s = sign[i]
            if s is None:
                y_vars.append([i, 1])
                y_vars.append([i, -1])
            elif s > 0:
                y_vars.append([i, 1])
            elif s < 0:
                y_vars.append([i, -1])
    return y_vars

def build_augmented_matrix(B, y_vars, A=None):
    m_B, n = B.shape
    n_y_vars = len(y_vars)
    
    M = np.concatenate((B, np.zeros((m_B, n_y_vars), dtype=int)), axis=1)
    for j in range(n_y_vars):
        i = y_vars[j][0]
        M[i][n+j] = -1*y_vars[j][1]
    
    if A is not None:
        m_A = A.shape[0]
        A_aug = np.concatenate((A, np.zeros((m_A, n_y_vars), dtype=int)), axis=1)
        M = np.concatenate((A_aug, M), axis=0)
        
    return M

def poly_circs(B_ineq, A_eq=None, sign=None):
    A = A_eq
    B = B_ineq
    m_B,n = B.shape
    m_a = 0
    if A is not None:
        m_a = A.shape[0]
    
    y_vars = determine_y_vars(m_B, sign=sign)
    n_y_vars = len(y_vars)
    
    #build constraint matrix M for conic polyhedral model Mr >= 0,
    #where first column of M is the r.h.s. vector 0.
    M = build_augmented_matrix(B, y_vars, A=A)
    M = np.concatenate((M, -1*M))
    I_y = np.concatenate((np.zeros((n_y_vars, n)), np.eye(n_y_vars)), axis=1).astype(int)
    M = np.concatenate((M, I_y))
    M = np.concatenate((np.zeros((2*m_a + 2*m_B + n_y_vars, 1), dtype=int), M), axis=1)
    
    #use cdd to enumerate extreme rays of the cone
    mat = cdd.Matrix(M, number_type='fraction')
    # mat = cdd.Matrix(M, number_type = 'float')
    mat.rep_type = cdd.RepType.INEQUALITY
    poly = cdd.Polyhedron(mat)
    rays = np.array(poly.get_generators())
    # print(rays)

    # obtain circuits from extreme rays
    circuits = []
    num_rays = rays.shape[0]
    for i in range(num_rays):
        g = rays[i,1:n+1]
        if not np.array_equal(g, np.zeros(n)):
            g = g*sympy.lcm([g[j].denominator for j in range(n) if g[j] != 0]) #normalize
            circuits.append(g)  
    return circuits