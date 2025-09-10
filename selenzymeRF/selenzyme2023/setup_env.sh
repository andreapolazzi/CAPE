# scripts/setup_env.sh
#!/usr/bin/env bash
set -euo pipefail

ENV_NAME=${1:-selenzymeRF}
YML=${2:-env/environment.yml}

# Create or update the env
if conda env list | grep -qE "^\s*${ENV_NAME}\s"; then
  conda env update -n "$ENV_NAME" -f "$YML"
else
  conda env create -n "$ENV_NAME" -f "$YML"
fi

eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

# Install PyTorch CPU wheels in a platform-aware way
OS="$(uname -s || echo Unknown)"
if [[ "$OS" == "Linux" ]]; then
  # Linux: use the CPU wheel index (fast, no CUDA)
  pip install --index-url https://download.pytorch.org/whl/cpu \
    torch==2.6.0
elif [[ "$OS" == "Darwin" ]]; then
  # macOS: regular PyPI works for CPU wheels
  pip install torch==2.6.0
else
  echo "Detected $OS. If this is Windows, please use WSL2 + Ubuntu and re-run."
  exit 1
fi

# Smoke-test the imports you actually use
python - <<'PY'
mods = ["rdkit","Bio","rxnmapper","rxn","pandas","numpy","scipy","sklearn","networkx","requests","tqdm","torch"]
missing=[]
for m in mods:
    try: __import__(m)
    except Exception as e: missing.append((m,str(e)))
print("Missing:", missing if missing else "none")
PY

echo "[OK] Environment '$ENV_NAME' ready."