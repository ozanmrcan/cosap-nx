# COSAP-NX Setup Guide

This guide covers the complete installation of COSAP-NX and its dependencies.

---

## Prerequisites

COSAP-NX requires:

| Dependency | Version | Purpose |
|------------|---------|---------|
| Python | >= 3.8 | Pipeline builder API |
| Nextflow | >= 23.04 | Workflow execution engine |
| Docker | any recent | Container runtime for bioinformatics tools |
| Java | >= 11 | Required by Nextflow |

---

## 1. Check Your Environment

Run these commands to verify what's already installed:

```bash
# Python
python3 --version

# Nextflow
nextflow -version

# Docker
docker --version

# Java (required for Nextflow)
java -version
```

---

## 2. Install Missing Dependencies

### 2.1 Java (if not installed)

Nextflow requires Java 11 or later.

```bash
# Ubuntu/Debian
sudo apt update
sudo apt install openjdk-17-jre-headless

# Verify
java -version
```

### 2.2 Nextflow (if not installed)

```bash
# Download and install
curl -s https://get.nextflow.io | bash

# Move to a location in your PATH
sudo mv nextflow /usr/local/bin/

# Verify
nextflow -version
```

### 2.3 Docker (if not installed)

```bash
# Ubuntu/Debian
sudo apt update
sudo apt install docker.io

# Add your user to the docker group (to run without sudo)
sudo usermod -aG docker $USER

# IMPORTANT: Log out and log back in for group changes to take effect
# Or run: newgrp docker

# Verify
docker run hello-world
```

---

## 3. Set Up Python Virtual Environment

### Why a Virtual Environment?

Modern Linux distributions protect the system Python from user packages (PEP 668). A virtual environment creates an **isolated Python installation** just for COSAP-NX, keeping your system clean and avoiding conflicts.

### 3.1 Create the Virtual Environment

```bash
# Navigate to the project
cd ~/grad_project/cosap-nx

# Create a virtual environment named 'venv'
python3 -m venv venv
```

This creates a `venv/` folder containing:
- A copy of the Python interpreter
- Its own `pip` for installing packages
- Isolated `site-packages` directory

### 3.2 Activate the Virtual Environment

```bash
# Activate (you'll see "(venv)" in your prompt)
source venv/bin/activate
```

Your prompt will change to show the active environment:
```
(venv) ozan@ozan-lenovo:~/grad_project/cosap-nx$
```

**Important:** You must activate the venv every time you open a new terminal to work on COSAP-NX.

### 3.3 Deactivate (when done working)

```bash
deactivate
```

---

## 4. Install COSAP-NX

With the virtual environment **activated**:

```bash
# Make sure you're in the cosap-nx directory
cd ~/grad_project/cosap-nx

# Activate venv (if not already active)
source venv/bin/activate

# Install COSAP-NX in development mode
pip install -e .
```

The `-e` flag means "editable" - changes you make to the code are immediately reflected without reinstalling.

### Verify Installation

```bash
# Check that the package is installed
pip list | grep cosap

# Try importing it
python -c "from cosap_nx import BamReader, VariantCaller, Pipeline; print('OK')"
```

---

## 5. Verify Complete Setup

Run this verification script:

```bash
# Make sure venv is activated
source venv/bin/activate

# Run verification
python -c "
import shutil
import subprocess

print('=== COSAP-NX Environment Check ===')
print()

# Python
import sys
print(f'Python: {sys.version}')

# COSAP-NX
try:
    import cosap_nx
    print(f'COSAP-NX: {cosap_nx.__version__}')
except ImportError as e:
    print(f'COSAP-NX: NOT INSTALLED - {e}')

# Nextflow
nf = shutil.which('nextflow')
if nf:
    result = subprocess.run(['nextflow', '-version'], capture_output=True, text=True)
    version = result.stdout.split('\\n')[0] if result.stdout else 'unknown'
    print(f'Nextflow: {version}')
else:
    print('Nextflow: NOT FOUND')

# Docker
docker = shutil.which('docker')
if docker:
    result = subprocess.run(['docker', '--version'], capture_output=True, text=True)
    print(f'Docker: {result.stdout.strip()}')
else:
    print('Docker: NOT FOUND')

print()
print('=== Ready to run pipelines! ===')
"
```

---

## 6. Quick Reference

### Daily Workflow

```bash
# 1. Navigate to project
cd ~/grad_project/cosap-nx

# 2. Activate virtual environment
source venv/bin/activate

# 3. Work on COSAP-NX...
python examples/germline_deepvariant.py

# 4. When done
deactivate
```

### Common Commands

```bash
# Activate venv
source venv/bin/activate

# Install/update COSAP-NX after code changes
pip install -e .

# Run a pipeline
python your_pipeline.py

# Run Nextflow directly (for testing)
nextflow run nf/main.nf -params-file params.json -profile docker

# Deactivate venv
deactivate
```

---

## Troubleshooting

### "externally-managed-environment" error

You forgot to activate the virtual environment:
```bash
source venv/bin/activate
```

### "nextflow: command not found"

Nextflow isn't in your PATH. Either:
- Install it: `curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/`
- Or add its location to PATH

### "permission denied" when running Docker

Your user isn't in the docker group:
```bash
sudo usermod -aG docker $USER
# Then log out and log back in
```

### Nextflow can't pull Docker images

Check Docker is running:
```bash
sudo systemctl start docker
docker run hello-world
```

---

## Next Steps

Once setup is complete, see:
- [User Guide](USER_GUIDE.md) - How to use COSAP-NX
- [API Reference](API_REFERENCE.md) - Detailed API documentation
