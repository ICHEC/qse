import numpy as np
import qat.qpus

import qse


def qse_test():
    qbits = qse.Qbits(positions=np.ones((4, 3)))
