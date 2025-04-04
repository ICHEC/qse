import qse 
import pulser
import numpy as np

qbits = qse.Qbits(positions=np.ones((4, 3)))


print("qse", qse.__version__)
print("np", np.__version__)
print("pulser", pulser.__version__)
print("All ok!")
