cimport numpy as np
import numpy as np
cimport cython

cdef extern from "math.h":
    double sqrt(double m)

@cython.boundscheck(False)
cpdef int compute_num_of_vecs(np.ndarray[np.double_t, ndim=2] prims,
                              np.ndarray[np.double_t, ndim=2] basis,
                              int nghbp, double dmax, int nq):
    """Computes the number of vectors for a given dmax value."""
    ncrq = 0
    #
    cdef int jq, l, m, n
    cdef double dx, sx, sy, sz
    for jq in range(nq):
        qmqpx=basis[jq,0] - basis[0,0]
        qmqpy=basis[jq,1] - basis[0,1]
        qmqpz=basis[jq,2] - basis[0,2]
        for l in range(-nghbp,nghbp+1):
            for m in range(-nghbp,nghbp+1):
                for n in range(-nghbp,nghbp+1):
                    sx = qmqpx+l*prims[0,0]+m*prims[1,0]+n*prims[2,0]
                    sy = qmqpy+l*prims[0,1]+m*prims[1,1]+n*prims[2,1]
                    sz = qmqpz+l*prims[0,2]+m*prims[1,2]+n*prims[2,2]
                    dx = sqrt(sx*sx+sy*sy+sz*sz)
                    if dx <= dmax:
                        ncrq += 1
    return ncrq
