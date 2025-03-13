#!/bin/env python3

import numpy as np
import qse
from qse.utils import squarelattice, triangularlattice, hexagonlattice, kagomelattice

# # square lattice
# square_reg = squarelattice(1, 6, 6)
# square_reg.draw(square_reg, radius=10)
# plt.show()

# # triangular lattice
# triangular_reg = triangularlattice(1,6,5)
# triangular_reg.draw(triangular_reg, radius=1.1)
# plt.show()

# # hexagonal lattice
# hexagonal_reg = hexagonlattice(1,5,5)
# hexagonal_reg.draw(hexagonal_reg, radius=1.1)
# plt.show()

# # kagome lattice
# kagome_reg = kagomelattice(1,5,5)
# kagome_reg.draw(kagome_reg, radius=3)
# plt.show()

# # construct lattice via lattice vectors
# N1, N2, N3 = 6, 6, 6
# a = 1 # lattice spacing
# # simple cubic unit cell
# sc_cell = qse.cell.Cell.fromcellpar([a, a, a, 90, 90, 90])

# # qbits object with one qubit at origin
# q_unit = qse.Qbits(positions=np.zeros((1,3)))

# # set q_unit's cell to simple cubic
# q_unit.cell = sc_cell

# # now repeat the pattern to get lattice of desired size
# q_sc = q_unit.repeat((N1, N2, N3))

# print(q_sc)

# q_sc.draw()
# plt.show()
# #fig = plt.figure()
# #ax = fig.add_subplot(projection='3d')
# #ax.scatter(q_sc.positions[:, 0], q_sc.positions[:, 1], q_sc.positions[:, 2], marker='.')
# #plt.show()
