import datetime
import sys

import numpy as np

import qse

_, calculator = sys.argv

# Define the lattice
repeats = 4
qbits = qse.lattices.chain(4.0, repeats)

duration = 400
omega0 = 10.01
delta0 = 0.12

amplitude = qse.Signal(np.ones(6) * omega0, duration)
detuning = qse.Signal(np.ones(6) * delta0, duration)

if calculator == "pulser":
    print("Running pulser calculation")
    calc = qse.calc.Pulser(
        amplitude=amplitude,
        detuning=detuning,
        qbits=qbits,
        label="test_run",
    )
    # Compute
    calc.build_sequence()

elif calculator == "myqlm":
    print("Running myqlm calculation")
    calc = qse.calc.Myqlm(
        amplitude=amplitude,
        detuning=detuning,
        qbits=qbits,
        label="test_run",
    )

calc.calculate()
spins = calc.get_spins()

results = {
    "calculator": calculator,
    "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    "qse": qse.__version__,
    "spins": spins,
    "state": calc.statevector,
}
np.save(f"results_{calculator}.npy", results)

for key, val in results.items():
    print(f"\n{key}:")
    print(val)
