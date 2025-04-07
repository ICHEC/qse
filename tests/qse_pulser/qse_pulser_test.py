import numpy as np
import pulser

import qse


def myqlm_test():
    qbits = qse.Qbits(positions=np.ones((4, 3)))
    print("pulser", pulser.__version__)
