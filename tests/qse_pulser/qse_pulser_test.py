import numpy as np
import pulser

import qse


def test_pulser():
    qbits = qse.Qbits(positions=np.ones((4, 3)))
    print("pulser", pulser.__version__)
    assert True
