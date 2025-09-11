# scripts/setup_env.ps1
# PowerShell translation of scripts/setup_env.sh
# Usage: .\setup_env.ps1 [-EnvName selenzymeRF] [-Yml env\environment.yml]
# Notes: Run this from an "Anaconda PowerShell Prompt" or a shell where Conda is initialized.

param(
    [string]$EnvName = "selenzymeRF",
    [string]$Yml = "env/environment.yml"
)

# Stop on first error (like `set -euo pipefail`)
$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Require-Command($name) {
    if (-not (Get-Command $name -ErrorAction SilentlyContinue)) {
        throw "Required command '$name' not found. Please install it or ensure it's on PATH."
    }
}

# Ensure conda exists
Require-Command conda

# Create or update the env (similar to `conda env list | grep ...`)
$envList = conda env list | Out-String
if ($envList -match "^\s*$([regex]::Escape($EnvName))\s") {
    Write-Host "Updating existing conda env '$EnvName' from $Yml ..."
    conda env update -n $EnvName -f $Yml
} else {
    Write-Host "Creating conda env '$EnvName' from $Yml ..."
    conda env create -n $EnvName -f $Yml
}

# Activate env
Write-Host "Activating env '$EnvName' ..."
conda activate $EnvName

# Install PyTorch CPU wheels in a platform-aware way
$platform = if ($env:OS) { $env:OS } else { [System.Environment]::OSVersion.Platform.ToString() }
Write-Host "Detected platform: $platform"

if ($platform -match "Windows") {
    # Windows: prefer the official CPU index
    pip install --index-url https://download.pytorch.org/whl/cpu torch==2.6.0
}
elseif ($platform -match "Linux" -or $platform -match "Unix") {
    # Linux: official CPU index
    pip install --index-url https://download.pytorch.org/whl/cpu torch==2.6.0
}
elseif ($platform -match "Darwin" -or $platform -match "Mac") {
    # macOS: CPU wheels available on PyPI
    pip install torch==2.6.0
}
else {
    Write-Warning "Unrecognized OS. Attempting CPU wheel install from the CPU index."
    pip install --index-url https://download.pytorch.org/whl/cpu torch==2.6.0
}

# Smoke-test the imports actually used
$py = @'
mods = ["rdkit","Bio","rxnmapper","rxn","pandas","numpy","scipy","sklearn","networkx","requests","tqdm","torch"]
missing=[]
for m in mods:
    try:
        __import__(m)
    except Exception as e:
        missing.append((m, str(e)))
print("Missing:", missing if missing else "none")
'@

# Execute the Python smoke test
$tmp = [System.IO.Path]::GetTempFileName()
[System.IO.File]::WriteAllText($tmp, $py)
try {
    python $tmp
} finally {
    Remove-Item $tmp -ErrorAction SilentlyContinue
}

Write-Host "[OK] Environment '$EnvName' ready."
