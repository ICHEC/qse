import numpy as np
import qse
from qse.calc.myqlm import MyQLM


def test_qse():
    qbits = qse.Qbits(positions=np.ones((4, 3)))
    assert True
