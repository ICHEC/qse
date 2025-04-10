It is recommmended to install within a new Conda environment.
To create a new conda environment use
```bash
conda create -n <my-env> python=<python_version>
```
where we recommend 3.10, 3.11 or 3.12 for the version of Python.
The environment can then be activated with
```bash
conda activate <my-env>
```

Make sure you're in the directory where the `pyproject.toml` file is situated and install the package and its dependencies using `pip`
```bash
pip install .
```
If you want to be able to edit the code and have the changes immediately take place without requiring a new installation, use the editable flag
```bash
pip install -e .
```
You can verify the package and dependencies were installed correctly with
```bash
pip list
```
if you installed `qse` in editable mode, `pip list` should show it having an editable project location.