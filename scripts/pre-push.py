import pathlib
import subprocess
import sys

# --- CONFIGURATION ---
SRC = ["qse", "scripts"]  # Only check your new Python code
PYPROJECT_PATH = pathlib.Path("pyproject.toml")


def run_check(name, command):
    """Runs a shell command. If it fails, the whole script exits."""
    print(f"--- Running {name} ---")
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        print(f"❌ {name} failed. Push aborted.")
        sys.exit(1)
    print(f"✅ {name} passed.")


def main():

    # Step 1: Ruff Check
    run_check("Ruff Check:", f"uv run ruff check {' '.join(SRC)}")

    # Step 2: Ruff Format
    run_check("Ruff Format:", f"uv run ruff format --check {' '.join(SRC)}")

    # Step 3: Run your tests (uncomment when ready)
    # run_check("Unit Tests", "python tests/simple-test.py")

if __name__ == "__main__":
    main()