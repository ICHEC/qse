#!/usr/bin/env python3
import re
from pathlib import Path

from tomlkit import dumps, parse

file_path = Path("pyproject.toml")

# Read and parse pyproject.toml
text = file_path.read_text(encoding="utf-8")
doc = parse(text)

# Get and increment version
version_str = doc["project"]["version"]
major, minor, patch = map(int, version_str.split("."))
new_version = f"{major}.{minor}.{patch + 1}"

# Update version in document
doc["project"]["version"] = new_version

# Save changes
file_path.write_text(dumps(doc), encoding="utf-8")

print(f"Bumped version: {version_str} → {new_version}")
