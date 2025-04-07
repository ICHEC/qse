import numpy as np
import qat.qpus

import qse
from qse.calc.pulser import MyQLM


def test_qse():
    qbits = qse.Qbits(positions=np.ones((4, 3)))
    assert True
