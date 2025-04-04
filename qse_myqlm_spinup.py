import qse 
import qat.qpus
import qlmaas.qpus
import numpy as np

qbits = qse.Qbits(positions=np.ones((4, 3)))


print("qse", qse.__version__)
print("np", np.__version__)
print("qat.qpus", qat.qpus.__version__)
print("qlmaas.qpus", qlmaas.qpus.__version__)
print("All ok!")
