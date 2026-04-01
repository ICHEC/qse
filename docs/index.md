---
jupytext:
    formats: md:myst
    text_representation:
        extension: .md
        format_name: myst
kernelspec:
    display_name: Python 3
    language: python
    name: python3
mystnb:
    render_markdown_format: myst
---

# QSE

The Quantum Simulation Environment (QSE) package is intended to provide a flexible,
modular way to frame a quantum simulation problem involving position-dependent quantum degrees of freedom.
Primarily we look at a collection of qubits at given set of coordinates. These can be arbitrary coordinates, or defined on a lattice.

QSE's design is adapted from Atomic Simulation Environment (ASE) to suit the needs
for an abstract representation for

1. **defining quantum computing systems**
2. **computing operations/simulations**

in a vendor and backend agnostic way. ASE's modular nature, and extensibility make it very useful for a similar quantum computing application.
Below is the visual organization of the components of QSE.

```{mermaid}
:config: {"layout": "tidy-tree"}

mindmap
   root((QSE))
      Qbits
         Qbit
         Cell
      Calculator
         Pulser
         MyQLM
         Qutip
      Utils
         Magnetic
         Signals
      lattices
         2D
         3D
      Visualise
```


```{note}

This project is under active development.
```



## Qbits

`Qbit` is the smallest class that represents just one qubits. `Qbits` is the primary class that represents a collection of qubits.
It can be instantiated by providing a list of coordinates, or as an empty class.
See the [Qbits examples](https://github/ICHEC/qse/docs/build/html/tutorials/creating_and_manipulating_qbits.html) for more details.



```{code-cell}
import qse
qsqr = qse.lattices.square(
   lattice_spacing=2.0,
   repeats_x=4, repeats_y=4)
qsqr.draw(radius=5.0)
```

## Calculator

Calculator

## Contributing

See the [contributing page](https://github.com/ICHEC/qse/blob/main/CONTRIBUTIONS.md).

