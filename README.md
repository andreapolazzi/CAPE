# iGEM Project - Setup

## Prerequisiti
- Installa **Miniconda**: https://docs.conda.io/en/latest/miniconda.html

## Setup e run (PUNTI 1 e 2 SOLO LA PRIMA VOLTA)
```bash
# 1. Clona il progetto
git clone [URL_REPO]
cd iGEM

# 2. Crea ambiente conda
conda env create -f environment.yml

# 3. Attiva ambiente
conda activate igem

# 4. Vai nella cartella Django
cd django_tool

# 5. Avvia il server
python manage.py runserver
```


## Per uscire
```bash
conda deactivate
```