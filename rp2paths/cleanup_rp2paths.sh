#!/bin/bash

echo "⏳ Disattivo eventuali ambienti Conda attivi..."
conda deactivate || true

echo "🧹 Rimuovo ambiente rp2paths_env (se esiste)..."
conda remove -n rp2paths_env --all -y || true

echo "🧹 Elimino eventuale cartella Micromamba..."
rm -rf ~/micromamba

echo "🧹 Ripulisco modifiche al file .zshrc (alias e variabili micromamba)..."
sed -i '' '/MAMBA_ROOT_PREFIX/d' ~/.zshrc
sed -i '' '/micromamba\/bin/d' ~/.zshrc

echo "✅ Pulizia completata."
echo "ℹ️ Ricorda: per rendere effettive le modifiche, riavvia il terminale o esegui: source ~/.zshrc"
