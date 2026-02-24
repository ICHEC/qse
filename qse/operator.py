import qutip as qp


class Operator:
    """
    Represents an operator in a quantum system.

    Parameters
    ----------
    operator: str | list[str]
        The type of qubit operator.
        Currently only "X", "Y", "Z" are supported. "N" the number operator is
        supported for QuTiP. If a list, must be equal in length to the size of
        the qubits tuple.
    qubits: int | list[int]
        A single integer or list of integers representing the qubits
        the operator acts on.
    nqubits: int
        The total number of qubits in the system.
    coef: float
        The coefficient associated with the term.
        Defaults to 1.

    Examples
    --------
    Create operator "XII"

    >>> qse.Operator("X", 0, 3)
    ... 1.00 X0

    Create operator "IIYIZ"

    >>> qse.Operator(["Y", "Z"], [2, 4], 5)
    ... 1.00 Y2 Z4
    """

    def __init__(self, operator, qubits, nqubits, coef=1.0):
        if isinstance(qubits, int):
            qubits = [qubits]
        if isinstance(operator, str):
            operator = [operator] * len(qubits)

        _check_operator(operator)

        if len(qubits) != len(operator):
            raise Exception(
                "The number of passed qubits must equal the number of passed operators."
            )

        self.operator = operator
        self.qubits = qubits
        self.nqubits = nqubits
        self.coef = coef

    def to_str(self):
        """
        Generates a string representation of the operator.

        The string is constructed as a tensor product of identity (I) and operator,
        where the operator is placed at the positions specified by `qubits`.

        Returns
        -------
        str
            A string of length `nqubits`, with "I" at all positions except for the
            qubits in `qubits`, which are replaced by the operator.
        """
        op = ["I"] * self.nqubits
        for qi, op_str in zip(self.qubits, self.operator):
            op[qi] = op_str
        return "".join(op)

    def to_qutip(self):
        """
        Generates a QuTiP representation of the operator.

        Returns
        -------
        qutip.Qobj
            The QuTiP operator.
        """
        op = [qp.qeye(2)] * self.nqubits
        for qi, op_str in zip(self.qubits, self.operator):
            op[qi] = _qutip_converter(op_str)
        return self.coef * qp.tensor(op)

    def __repr__(self):
        return f"{self.coef:.2f} " + " ".join(
            [f"{op}{q}" for op, q in zip(self.operator, self.qubits)]
        )


def _check_operator(op_list):
    for op in op_list:
        if op not in ["X", "Y", "Z", "N"]:
            raise Exception("An operator must be one of 'X', 'Y', 'Z' or 'N'.")


def _qutip_converter(op):
    if op == "X":
        return qp.sigmax()

    if op == "Y":
        return qp.sigmay()

    if op == "Z":
        return qp.sigmaz()

    if op == "N":
        return 0.5 * (1 - qp.sigmaz())
