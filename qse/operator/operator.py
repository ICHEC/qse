import qutip as qp


class Operator:
    """
    Represents an operator in a quantum system.

    Parameters
    ----------
    operator: str | list[str]
        The type of qubit operator.
        Currently only a string or list containing "X", "Y", "Z", "N" are supported.
        "N" the number operator is defined as 0.5*(1-Z).
        If a list is passed, it must be equal in length to the size of
        the qubits list.
    qubits: int | list[int]
        A single integer or list of integers representing the qubits
        the operator acts on.
    nqbits: int
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

    def __init__(self, operator, qubits, nqbits, coef=1.0):
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
        self.nqbits = nqbits
        self.coef = coef

    def to_str(self):
        """
        Generates a string representation of the operator.

        The string is constructed as a tensor product of identity (I) and operator,
        where the operator is placed at the positions specified by `qubits`.

        Returns
        -------
        str
            A string of length `nqbits`, with "I" at all positions except for the
            qubits in `qubits`, which are replaced by the associated operator.
        """
        op = ["I"] * self.nqbits
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
        op = [qp.qeye(2)] * self.nqbits
        for qi, op_str in zip(self.qubits, self.operator):
            op[qi] = _qutip_converter(op_str)
        return self.coef * qp.tensor(op)

    def extend(self, other):
        """
        Extend Operator object by appending terms from other SumOfOperators
        object or Operator object.

        Parameters
        ----------
        other: Operators
            The operators to be added to the current object.
        """
        if other.nqbits != self.nqbits:
            raise Exception("other must have the same number of qubits.")

        if isinstance(other, Operator):
            self.operator_list += [other]
        elif isinstance(other, SumOfOperators):
            self.operator_list += other.operator_list
        else:
            raise Exception("other must be SumOfOperators or Operator.")

    def __add__(self, other):
        op = self.copy()
        op += other
        return op

    def __iadd__(self, other):
        self.extend(other)
        return self

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
