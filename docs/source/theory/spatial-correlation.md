# Computing spatial correlation in magnetic models

Imagine a lattice of say $N$ atoms, in $L_1\times L_2\times L_3$ geometry of repititions,
with ${\bf A_1, A_2, A_3}$ as lattice vectors. A general point on this lattice ie given by

$$
{\bf X(n_1, n_2, n_3)} = n_1 {\bf A_1} + n_2 {\bf A_2} + n_3 {\bf A_3}
$$

Where $n_1, n_2, n_3$ are integers, and for the finite size chosen above, $n_i = 0,1,\dots L_i - 1$.
We shorten the labels as $i\equiv (n_1, n_2, n_3)$, and express lattice positions as ${\bf X_i}$.

## Classical Spins

Now, imagine a lattice where each of the site has a classical spin, represented by a unit vector ${\bf S_i}$. If the energy of this system is given by a classical Heisenberg model, then we have

$$
H = \sum_{ij} J_{ij}~{\bf S_i}\cdot {\bf S_j}
$$

where $J_{ij}$ is exchange between the spins, often restricted to nearest neighbours or next nearest neighbours.

Depending on the exchange, the ground state, or lowest energy state of the system will have some specific configuration of spin vectors. For example, if $J_{ij} < 0$, then independent of the scale of exchanges, the energy is minimized in ferromagnetic state, i.e., when all the spins are pointed in the same direction ${\bf S_i} = (0, 0, 1)$.


When spins are in states other than ferromagnetic state, one needs some quantity for better comprehension of the nature of the spin orientation. There are several commonly known periodic states that can be visualised spatially, such as Neel type Antiferromagnetic states. Howevever, more complex periodic patterns, or somewhat disordered patterns are hard to make sense of visually. For this we compute the following:

### Spin-Correlation

One computes $S_{ij} = \langle {\bf S_i} \cdot {S_j}\rangle$ for all $i, j$. For classical models, this could be for a single spin configuration, or at finite temperature, averaged of thermally generated ensemble of configurations. Then we compute Fourier transform of this: 

$$
S({\bf q}) = \frac {1}{N^2} \sum_{ij} S_{ij}~e^{\left(i {\bf q}\cdot ({\bf X_i} - {\bf X_j})\right)}
$$

This is defined as the **Structure Factor**.

### Properties

To be added

### Examples

To be added

## Quantum Spins

Let's consider case of quantum spins, that is, say at each site we have a spin-$\frac 12$ object instead of classical spins. In that case we have to know the many-body quantum state, say $|\Psi \rangle$ in which the system exists, and the site correlation we need to compute is the expectation value of the spin operators:

$$
\langle\Psi| {\bf S_i}\cdot {\bf S_j}|\Psi\rangle
$$

To compute the expectation value $\langle \Psi | \mathbf{S}_i \cdot \mathbf{S}_j | \Psi \rangle$ for a system of $N$ spin-$\frac 12$ particles, where $\mathbf{S}_i = (S_i^x, S_i^y, S_i^z)$ are the spin operators at site $i$, we derive the expression by expanding the dot product and simplifying using the  properties of Pauli matrices. The spin operators are defined as $S_i^\alpha = \frac{\hbar}{2} \sigma_i^\alpha$ for $\alpha = x, y, z$, where $\sigma_i^\alpha$ are the Pauli matrices. Setting $\hbar = 1$ for simplicity, we have $S_i^\alpha = \frac{1}{2} \sigma_i^\alpha$.

The operator $\mathbf{S}_i \cdot \mathbf{S}_j$ is:

$$
\mathbf{S}_i \cdot \mathbf{S}_j = S_i^x S_j^x + S_i^y S_j^y + S_i^z S_j^z.
$$

Substituting the expressions in terms of Pauli matrices:

$$
\mathbf{S}_i \cdot \mathbf{S}_j = \frac{1}{4} \left( \sigma_i^x \sigma_j^x + \sigma_i^y \sigma_j^y + \sigma_i^z \sigma_j^z \right).
$$

We now use the identity for the sum of the \(x\) and \(y\) components:

$$
\sigma_i^x \sigma_j^x + \sigma_i^y \sigma_j^y = \sigma_i^+ \sigma_j^- + \sigma_i^- \sigma_j^+,
$$

where $\sigma_i^+ = \sigma_i^x + i\sigma_i^y$ and $\sigma_i^- = \sigma_i^x - i\sigma_i^y$ are the raising and lowering operators. This identity holds for all $i$ and $j$, but we must handle the case $i = j$ separately due to the properties of the operators.

### Case 1: $i = j$
When $i = j$, the expectation value simplifies to:

$$
\mathbf{S}_i \cdot \mathbf{S}_i = (S_i^x)^2 + (S_i^y)^2 + (S_i^z)^2.
$$

For a spin-$\frac 12$ particle, the eigenvalues of $\mathbf{S}_i^2$ are $s(s+1) = \frac{3}{4}$ (since $s = \frac{1}{2}$), and the operator is proportional to the identity:

$$
\mathbf{S}_i \cdot \mathbf{S}_i = \frac{3}{4} \mathbb{I}.
$$

Thus, the expectation value is:

$$
\langle \Psi | \mathbf{S}_i \cdot \mathbf{S}_i | \Psi \rangle = \frac{3}{4} \langle \Psi | \mathbb{I} | \Psi \rangle = \frac{3}{4}.
$$

### Case 2: $i \neq j$
For $i \neq j$, we use the expanded form:

$$
\mathbf{S}_i \cdot \mathbf{S}_j = \frac{1}{4} \left( \sigma_i^x \sigma_j^x + \sigma_i^y \sigma_j^y + \sigma_i^z \sigma_j^z \right) = \frac{1}{4} \left( \sigma_i^+ \sigma_j^- + \sigma_i^- \sigma_j^+ + \sigma_i^z \sigma_j^z \right).
$$

Thus, the expectation value is:

$$
\langle \Psi | \mathbf{S}_i \cdot \mathbf{S}_j | \Psi \rangle = \frac{1}{4} \left[ \langle \Psi | \sigma_i^+ \sigma_j^- | \Psi \rangle + \langle \Psi | \sigma_i^- \sigma_j^+ | \Psi \rangle + \langle \Psi | \sigma_i^z \sigma_j^z | \Psi \rangle \right].
$$

### Final Expression
Combining both cases, the expectation value is:

$$
\langle \Psi | \mathbf{S}_i \cdot \mathbf{S}_j | \Psi \rangle = 
\begin{cases} 
\frac{3}{4} & \text{if } i = j, \\
\frac{1}{4} \left( \langle \sigma_i^+ \sigma_j^- \rangle + \langle \sigma_i^- \sigma_j^+ \rangle + \langle \sigma_i^z \sigma_j^z \rangle \right) & \text{if } i \neq j,
\end{cases}
$$

where $\langle \cdot \rangle$ denotes the expectation value in the state $|\Psi\rangle$.

### Summary
- For $i = j$, the result is $\frac{3}{4}$.
- For $i \neq j$, the expectation value requires computing three terms: the expectation values of the operators $\sigma_i^+ \sigma_j^-$, $\sigma_i^- \sigma_j^+$, and $\sigma_i^z \sigma_j^z$ in the state $|\Psi\rangle$.

This expression accounts for the quantum mechanical nature of the spin operators and is valid for any state $|\Psi\rangle$ in the Hilbert space. The operators $\sigma_i^\pm$ and $\sigma_i^z$ are defined in the basis of the $2^N$-dimensional Hilbert space, and their matrix elements can be evaluated once the state $|\Psi\rangle$ is specified.


## Derivation of the final expression

Now let's get a bit more in detail, and express the state $\Psi$ as linear combination of the product  state basis of the $N$-qubit or $N$-spin Hilbert space. Then $|\Psi\rangle = \sum_{\alpha}|\alpha\rangle$ where, $\alpha$ symbolically is $\alpha \equiv a_0a_1a_2 \dots a_{N-1}$, and $a_i$ are binary numbers.


To compute the expectation value $\langle \Psi | \mathbf{S}_i \cdot \mathbf{S}_j | \Psi \rangle$ for a state $|\Psi\rangle = \sum_a c_a |a\rangle$, where $|a\rangle$ are computational basis states (e.g., $|a\rangle = |00101\rangle$), we express the operator $\mathbf{S}_i \cdot \mathbf{S}_j$ in terms of Pauli matrices and evaluate its expectation using the coefficients $c_a$. The derivation below assumes the basis convention where:

- $|0\rangle$ corresponds to $\sigma^z = +1$ (spin up in $z$-direction).
- $|1\rangle$ corresponds to $\sigma^z = -1$ (spin down in $z$-direction).

### 1. **Operator Decomposition**
The spin operator $\mathbf{S}_i \cdot \mathbf{S}_j$ decomposes as:

$$
\mathbf{S}_i \cdot \mathbf{S}_j = \frac{1}{4} \left( \sigma_i^x \sigma_j^x + \sigma_i^y \sigma_j^y + \sigma_i^z \sigma_j^z \right).
$$

Using raising/lowering operators $(\sigma^\pm = \frac{1}{2}(\sigma^x \pm i\sigma^y)$):

$$
\sigma_i^x \sigma_j^x + \sigma_i^y \sigma_j^y = 2(\sigma_i^+ \sigma_j^- + \sigma_i^- \sigma_j^+),
$$

so the operator simplifies to:

$$
\mathbf{S}_i \cdot \mathbf{S}_j = \frac{1}{4} \left( 2\sigma_i^+ \sigma_j^- + 2\sigma_i^- \sigma_j^+ + \sigma_i^z \sigma_j^z \right) = \frac{1}{2} \sigma_i^+ \sigma_j^- + \frac{1}{2} \sigma_i^- \sigma_j^+ + \frac{1}{4} \sigma_i^z \sigma_j^z.
$$

### 2. **Expectation Value for $i = j$**
When $i = j$:

$$
\langle \mathbf{S}_i \cdot \mathbf{S}_i \rangle = \langle (\mathbf{S}_i)^2 \rangle = s(s+1) = \frac{3}{4}.
$$

### 3. **Expectation Value for $i \neq j$**
For $i \neq j$, the expectation value is:

$$
\langle \mathbf{S}_i \cdot \mathbf{S}_j \rangle = \frac{1}{2} \langle \sigma_i^+ \sigma_j^- \rangle + \frac{1}{2} \langle \sigma_i^- \sigma_j^+ \rangle + \frac{1}{4} \langle \sigma_i^z \sigma_j^z \rangle.
$$

We evaluate each term using the coefficients $c_a$:

#### (a) **Diagonal Term $(\langle \sigma_i^z \sigma_j^z \rangle$)**
This term is **diagonal** in the basis, so:

$$
\sigma_i^z \sigma_j^z |a\rangle = (1 - 2a_i)(1 - 2a_j) |a\rangle,
$$

where $a_i, a_j \in \{0, 1\}$ are the bits of $|a\rangle$ at positions $i$ and $j$. Thus:

$$
\langle \sigma_i^z \sigma_j^z \rangle = \sum_{a} |c_a|^2 (1 - 2a_i)(1 - 2a_j).
$$

#### (b) **Off-Diagonal Term 1 ($\langle \sigma_i^+ \sigma_j^- \rangle$)**
- $\sigma_j^-$ flips $a_j$ from 0 to 1 (if $a_j = 0$).
- $\sigma_i^+$ flips $a_i$ from 1 to 0 (if $a_i = 1$).
- Non-zero only when $a_j = 0$ and $a_i = 1$:

$$
\sigma_i^+ \sigma_j^- |a\rangle = |b\rangle, \quad \text{where } b = a \text{ with } a_j = 1 \text{ and } a_i = 0.
$$

The expectation is:

$$
\langle \sigma_i^+ \sigma_j^- \rangle = \sum_{\substack{a \\ a_i=1, a_j=0}} c_a^* c_b, \quad b = a \oplus (2^i + 2^j).
$$

#### (c) **Off-Diagonal Term 2 ($\langle \sigma_i^- \sigma_j^+ \rangle$)**
- $\sigma_j^+$ flips $a_j$ from 1 to 0 (if $a_j = 1$).
- $\sigma_i^-$ flips $a_i$ from 0 to 1 (if $a_i = 0$).
- Non-zero only when $a_j = 1$ and $a_i = 0$:

$$
\sigma_i^- \sigma_j^+ |a\rangle = |d\rangle, \quad \text{where } d = a \text{ with } a_j = 0 \text{ and } a_i = 1.
$$

The expectation is:

$$
\langle \sigma_i^- \sigma_j^+ \rangle = \sum_{\substack{a \\ a_i=0, a_j=1}} c_a^* c_d, \quad d = a \oplus (2^i + 2^j).
$$

### 4. **Final Expression for $i \neq j$**
Combining all terms, we have the following:

$$
\langle \mathbf{S}_i \cdot \mathbf{S}_j \rangle = \frac{1}{2} \left( \sum_{\substack{a \\ a_i=1, a_j=0}} c_a^* c_b \right) + \frac{1}{2} \left( \sum_{\substack{a \\ a_i=0, a_j=1}} c_a^* c_d \right) + \frac{1}{4} \left( \sum_{a} |c_a|^2 (1 - 2a_i)(1 - 2a_j) \right)
$$

where:
- $b = a \oplus (2^i + 2^j)$ (flip bits at $i$ and $j$),
- $d = a \oplus (2^i + 2^j)$ (same as $b$).

### Summary
- **For $i = j$**: $\langle \mathbf{S}_i \cdot \mathbf{S}_j \rangle = \frac{3}{4}$.
- **For $i \neq j$**: Use the boxed expression above, which depends on:
  1. Diagonal terms: Weighted by $|c_a|^2$ and spin eigenvalues.
  2. Off-diagonal terms: Correlations between basis states differing by spin flips at $i$ and $j$.

This expression is efficiently computable by iterating over all basis states $a$ and summing contributions.



### rough work

$$
|\alpha\rangle \equiv |a_0 a_1 a_2 ... a_N\rangle
 = |a_0\rangle \times |a_1\rangle \times |a_2\rangle \times ... |a_N\rangle
$$

$$
(p_0^0|+\rangle + p_1^0|-\rangle) \times
(p_0^1|+\rangle + p_1^1|-\rangle) \times
(p_0^2|+\rangle + p_1^2|-\rangle) \times
\dots
(p_0^N|+\rangle + p_1^N|-\rangle)
$$

$$
\langle S_{x,y,z}\rangle
$$

$$
\langle S_{x,y,z}S_{x,y,z}\rangle
$$
$$
|\psi\rangle = \sum_\alpha c_\alpha |\alpha \rangle
$$

