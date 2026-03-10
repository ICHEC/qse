from .operator import Operator


class Operators:
    """
    Represents a sum of operators in a quantum system.

    Parameters
    ----------
    operator_list: list[Operator] , optional
        The operators to be summed.
        Defaults to empty list.
    """

    def __init__(self, operator_list=None):
        if operator_list is None:
            operator_list = []
        else:
            if not isinstance(operator_list, list) and not isinstance(
                operator_list[0], Operator
            ):
                raise Exception("The operator_list must be a list of qse.Operators.")

            if not all([op.nqbits == operator_list[0].nqbits for op in operator_list]):
                raise Exception(
                    "All operators in operator_list must act on the same number"
                    "of qubits."
                )

        self.operator_list = operator_list

    def to_qutip(self):
        """
        Generates a QuTiP representation of the sum of operators.

        Returns
        -------
        qutip.Qobj
            The QuTiP operator.
        """
        operator = 0
        for op in self.operator_list:
            operator += op.to_qutip()
        return operator

    def to_qiskit(self):
        """
        Generates a qiskit representation of the sum of operators.

        Returns
        -------
        qiskit.quantum_info.SparsePauliOp
            The qiskit sparse pauli operator.
        """
        from qiskit.quantum_info import SparsePauliOp

        return SparsePauliOp(
            [i.to_str() for i in self.operator_list],
            [i.coef for i in self.operator_list],
        )

    def extend(self, other):
        """
        Extend Operators object by appending terms from other Operators
        object or Operator object.

        Parameters
        ----------
        other: Operators or Operator
            The operators to be added to the current object.
        """
        if other.nqbits != self.nqbits:
            raise Exception("other must have the same number of qubits.")

        if isinstance(other, Operator):
            self.operator_list += [other]
        elif isinstance(other, Operators):
            self.operator_list += other.operator_list
        else:
            raise Exception("other must be Operators or Operator.")

    def __add__(self, other):
        op = self.copy()
        op += other
        return op

    def __iadd__(self, other):
        self.extend(other)
        return self

    def __repr__(self):
        return (
            f"Number of qubits: {self.nqbits}"
            + f"\nNumber of terms: {self.size}\n\n"
            + "\n".join([op.__repr__() for op in self.operator_list])
        )

    def __getitem__(self, i):
        return self.operator_list[i]

    @property
    def nqbits(self):
        return self[0].nqbits

    @property
    def nterms(self):
        return len(self.operator_list)
