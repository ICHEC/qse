class InteractionTerm:
    """
    Represents an interaction term in a quantum system.

    This class encapsulates an interaction between two qubits, including the type of interaction,
    the qubits involved, the coefficient, and the total number of qubits in the system.

    Parameters
    ----------
    interaction: str 
        The type of interaction (e.g., "X", "Y", "Z").
    qubits: tuple
        A tuple of two integers representing the indices of the interacting qubits.
    coef: float
        The coefficient associated with the interaction term.
    nqubits: int
        The total number of qubits in the system.
    """
    def __init__(self, interaction, qubits, coef, nqubits):        
        self.interaction = interaction
        self.qubits = qubits
        self.coef = coef
        self.nqubits = nqubits

    def to_str(self):
        """
        Generates a string representation of the interaction term.

        The string is constructed as a tensor product of identity (I) and interaction operators,
        where the interaction operator is placed at the positions specified by `qubits`.

        Returns
        -------
        str
            A string of length `nqubits`, with "I" at all positions except for the qubits in `qubits`,
            which are replaced by the `interaction` operator.
        """
        op = ["I"] * self.nqubits
        for qi in self.qubits:
            op[qi] = self.interaction
        return "".join(op)
