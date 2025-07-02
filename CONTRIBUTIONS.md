# Contributing
We adhere to [PEP 8](https://peps.python.org/pep-0008/) style guidelines and use the following formatting tools:
- [Flake8](https://flake8.pycqa.org/en/latest/) (for style guide enforcement)
- [isort](https://pycqa.github.io/isort/) (for sorting imports)
- [black](https://github.com/psf/black) (for code formatting)

and we use [pytest](https://docs.pytest.org/en/stable/) for running unit tests.

To install these tools use
```bash
pip install flake8 black isort pytest
```
or (from inside the directory where the `pyproject.toml` is located)
```bash
pip install ".[dev]"
```

To run `black` on all python files within a directory, use
```bash
black .
```
Similarly to run `isort` on all python files within a directory, use
```bash
isort .
```
One can run `black` or `isort` on a single file by changing `.` to `<file_name>`.

All tests should be in the `tests` directory.
To run `pytest` use
```bash
pytest
```
this will search recursively for any files of the form `test_*.py` or `_test.py` (for information on how `pytest` searches for tests see [here](https://docs.pytest.org/en/7.1.x/explanation/goodpractices.html#conventions-for-python-test-discovery)). 
To run `pytest` on a single file use
```bash
pytest <filename>
```
>> [NOTE] For developing, if one clones the git repo, one should run `./setup.sh` script to setup the necessary git-hooks in the repo after cloning it.