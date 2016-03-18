from emto_input_creator import *

# Some test cases
#"""
# L11 CuPt
prims1 = np.array([[1.0,0.5,0.5],
                   [0.5,1.0,0.5],
                   [0.5,0.5,1.0]])
basis1 = np.array([[0.0,0.0,0.0],
                   [0.5,0.5,0.5]])
lattice1 = Lattice(prims1)
struct1 = Structure(lattice1, ["Cu","Pt"], basis1, coords_are_cartesian=False)
#"""
"""
# conventional fcc
prims1 = a*np.array([[1.0,0.0,0.0],
                     [0.0,1.0,0.0],
                     [0.0,0.0,1.0]])
basis1 = a*np.array([[0.0,0.0,0.0],
                     [0.5,0.5,0.0],
                     [0.5,0.0,0.5],
                     [0.0,0.5,0.5]])
"""
"""
# fcc
prims1 = a*np.array([[0.5,0.5,0.0],
                     [0.5,0.0,0.5],
                     [0.0,0.5,0.5]])
basis1 = np.array([[0,0,0]])
"""

"""
# bcc SSOS-1
prims1 = a*np.array([[0.0,2.0,3.0],
                     [-0.5,2.5,2.5],
                     [-0.5,-0.5,0.5]])

basis1 = np.array([[0.0,0.0,0.0],
                   [0.6,0.6,0.4],
                   [0.4,0.4,0.6],
                   [0.2,0.2,0.8],
                   [0.8,0.8,0.2]])
lattice1 = Lattice(prims1)
struct1 = Structure(lattice1, ["Fe","Cr","Mo","Ni","Cu"], basis1, coords_are_cartesian=False)
"""

"""
# bcc SSOS-2
prims1 = a*np.array([[-0.5,-1.5,-0.5],
                     [-0.5,1.5,0.5],
                     [0.5,0.5,-1.5]])

basis1 = np.array([[0.4,0.2,0.6],
                   [0.2,0.6,0.8],
                   [0.6,0.8,0.4],
                   [0.8,0.4,0.2],
                   [0.0,0.0,0.0]])
lattice1 = Lattice(prims1)
struct1 = Structure(lattice1, ["Fe","Cr","Mo","Ni","Cu"], basis1, coords_are_cartesian=False)
"""

"""
# bcc SSOS-3
prims1 = a*np.array([[-1.5,-0.5,-1.5],
                     [-1.5,-1.5,-0.5],
                     [-0.5,0.5,0.5]])

basis1 = np.array([[0.0,0.0,0.0],
                   [0.2,0.2,0.8],
                   [0.8,0.8,0.2],
                   [0.4,0.4,0.6],
                   [0.6,0.6,0.4]])
lattice1 = Lattice(prims1)
struct1 = Structure(lattice1, ["Fe","Cr","Mo","Ni","Cu"], basis1, coords_are_cartesian=False)
"""
"""
# fcc SSOS-1
prims1 = a*np.array([[-1.5,1.5,2.0],
                     [-1.5,2.0,1.5],
                     [-2.0,1.5,1.5]])

basis1 = np.array([[0.0,0.0,0.0],
                   [0.2,0.2,0.2],
                   [0.8,0.8,0.8],
                   [0.6,0.6,0.6],
                   [0.4,0.4,0.4]])
lattice1 = Lattice(prims1)
struct1 = Structure(lattice1, ["Fe","Cr","Mo","Ni","Cu"], basis1, coords_are_cartesian=False)
"""
"""
# fcc SSOS-2
prims1 = a*np.array([[-1.5,1.5,2.0],
                     [-1.5,2.0,1.5],
                     [-2.0,1.5,1.5]])

basis1 = np.array([[0.0,0.0,0.0],
                   [0.4,0.4,0.4],
                   [0.6,0.6,0.6],
                   [0.2,0.2,0.2],
                   [0.8,0.8,0.8]])
lattice1 = Lattice(prims1)
struct1 = Structure(lattice1, ["Fe","Cr","Mo","Ni","Cu"], basis1, coords_are_cartesian=False)
"""
"""
# fcc SSOS-3
prims1 = a*np.array([[0.5,-1.0,0.5],
                     [-0.5,1.0,0.5],
                     [-1.0,-0.5,-0.5]])

basis1 = np.array([[0.0,0.0,0.0],
                   [0.8,0.4,0.2],
                   [0.2,0.6,0.8],
                   [0.4,0.2,0.6],
                   [0.6,0.8,0.4]])
lattice1 = Lattice(prims1)
struct1 = Structure(lattice1, ["Fe","Cr","Mo","Ni","Cu"], basis1, coords_are_cartesian=False)
"""
