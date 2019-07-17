from __future__ import division
import numpy as np
cimport numpy as np
cimport cython


DTYPE = np.float
ctypedef np.float_t DTYPE_t

# Matrix multiplication by a diagonal is equivalent to multiplying
# the columns of A by the nonzero elements of the columns of D.

@cython.boundscheck(False)
@cython.wraparound(False)
def compute_D(np.ndarray[DTYPE_t, ndim=1] P):
    cdef int Px = P.shape[0]
    cdef np.ndarray[DTYPE_t, ndim=1] D = np.empty(Px, dtype=np.float)
    cdef int x

    for x in range(Px):
        D[x] = (1/Px)*(2*P[x])*(1 - P[x])
    return D

@cython.boundscheck(False)
@cython.wraparound(False)
def diag_matrix_mult(np.ndarray[DTYPE_t, ndim=2] A, np.ndarray[DTYPE_t, ndim=1] D):
    cdef int Ax = A.shape[0]
    cdef int Ay = A.shape[1]
    cdef int Dx = D.shape[0]
    cdef np.ndarray[DTYPE_t, ndim=2] C = np.empty([Ax, Ay], dtype=np.float)
    cdef int x, y

    for y in range(Ay):
        for x in range(Ax):
            C[x, y] = A[x, y]*D[x]

    return C

@cython.boundscheck(False)
@cython.wraparound(False)
def compute_G(np.ndarray[DTYPE_t, ndim=2] ZD, np.ndarray[DTYPE_t, ndim=2] Z):
    cdef int ZDx, ZDy, Zx, Zy
    cdef int i, j

    ZDx, ZDy = ZD.shape[0], ZD.shape[1]
    Zx, Zy = Z.shape[0], Z.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=2] ZDZT = np.zeros([ZDx, ZDx], dtype=np.float)

    i = 0
    j = 0
    for i in range(ZDx):
        for j in range(ZDx):
            ZDZT[i, j] = np.sum(ZD[i, :]*Z.T[:, j])
    return ZDZT


