# Removal of unused functions/classes/methods

The following has been removed:

| function/class/method                                 | purpose/source   | Action    |
|-------------------------------------------------------|------------------|-----------|
| qse/constraints.py [all]                              |                  |           |
| qse/calc/messages.py [all]                            |                  |           |
| qse.calc.calculator.compare_qbits                     | ASE              | Remove    |
| qse.calc.calculator.equal                             | ASE              | Remove    |
| qse.calc.calculator.Parameters                        |                  |           |
| qse.calc.calculator.Calculator.get_default_parameters | sounds redundant | Remove    |
| qse.calc.calculator.Calculator.reset                  | Can be used      | Repurpose |
| qse.calc.calculator.Calculator.read                   | Can be used      | Repurpose |
| qse.calc.calculator.Calculator.get_qbits              |                  |           |
| qse.calc.calculator.Calculator.read_qbits             |                  |           |
| qse.calc.calculator.Calculator.set                    |                  |           |
| qse.calc.calculator.Calculator.check_state            |                  |           |
| qse.calc.calculator.Calculator.get_energy             |                  |           |
| qse.calc.calculator.Calculator.get_property           |                  |           |
| qse.calc.calculator.Calculator.calculation_required   |                  |           |
| qse.calc.calculator.Calculator.export_properties      |                  |           |
| qse.calc.calculator.FileIOCalculator                  |                  |           |
| qse.qbits.celldisp                                    |                  |           |
| qse.qbits.Qbits.set_constraint                        |                  |           |
| qse.qbits.Qbits.constraint                            |                  |           |
| qse.qbits.Qbits.set_celldisp                          |                  |           |
| qse.qbits.Qbits.get_celldisp                          |                  |           |
| qse.qbits.Qbits.get_cell_lengths_and_angles           |                  |           |
| qse.qbits.Qbits.get_reciprocal_cell                   |                  |           |
| qse.qbits.Qbits.pbc                                   |                  |           |
| qse.qbits.Qbits.get_properties                        |                  |           |
| qse.qbits.Qbits.center_in_unit_cell                   |                  |           |
| qse.qbits.Qbits.get_dihedral                          |                  |           |
| qse.qbits.Qbits.get_dihedrals                         |                  |           |
| qse.qbits.Qbits.set_dihedral                          |                  |           |
| qse.qbits.Qbits.rotate_dihedral                       |                  |           |
| qse.qbits.Qbits.get_scaled_positions                  |                  |           |
| qse.qbits.Qbits.wrap                                  |                  |           |
