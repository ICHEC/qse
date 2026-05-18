import pathlib
import shutil
import subprocess
import sys

# --- CONFIGURATION ---
SRC = [
    ".",
]  # Only check your new Python code


def get_uv_prefix():
    """
    Determines if we should prefix commands with 'uv run'.
    Checks if 'uv' is installed and if a 'uv.lock' file exists.
    """
    has_uv = shutil.which("uv") is not None
    is_uv_project = pathlib.Path("uv.lock").exists()

    if has_uv and is_uv_project:
        return "uv run "
    return ""


def run_check(name, command):
    """
    Detects environment and runs the appropriate command.
    If it fails, the whole script exits.
    """
    prefix = get_uv_prefix()
    full_command = f"{prefix}{command}"

    print(f"---🔍 Running {name} ---")
    print(f"Executing: {full_command}")

    result = subprocess.run(command, shell=True)

    if result.returncode != 0:
        print(f"❌ {name} failed. Push aborted.")
        sys.exit(1)
    print(f"✅ {name} passed.")


def main():
    # Define the base commands without the 'uv run' part

    # Step 1: Ruff Check
    run_check("Ruff Check:", f"ruff check {' '.join(SRC)} --output-format grouped")

    # Step 2: Ruff Format
    run_check("Ruff Format:", f"ruff format --check {' '.join(SRC)}")

    # Step 3: Run your tests (uncomment when ready)
    # run_check("Unit Tests", "python tests/simple-test.py")
    print("\n🚀 All checks passed!")


if __name__ == "__main__":
    main()
