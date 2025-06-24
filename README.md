# Struttura della repo

```
cartella_privata/
|
├── iGEM
|   ├── codes/
|   ├── data_files/
|   └── results/
|
|
└── retrorules_rr02_rp2_flat_forward.csv
```

* Il file
  ```retrorules_rr02_rp2_flat_forward.csv```
  contiene le regole di reazione. Scaricabile da [RetroRules](https://retrorules.org/dl), sezione rr02 per RetroPath 2.0. Deve essere posto al di fuori della repo.
* iGEM è la repo vera e propria.
  * ```codes/``` contiene gli script
    * ```performance.py``` è usato per testare le risorse computazionali usate da retropath
    * ```results_parser.py``` è una prova di un parser dei risultati di retropath quando non produce scope, ma la stiamo sostituendo con una versione basata sui grafi. Al momento è inutilizzato
  * ```data_files/``` contiene i file di dati, come i vari sink e i source
    * ```biocyc_cmpd_rpcleaned.csv``` file con tutti i metaboliti di rhodococcus in formato utile per retropath
    * ```gem_central_rp.csv``` file con i metaboliti centrali di rhodococcus opacus in formato utile per retropath
    * ```log*.txt``` sono tutti i file con i log dei test che vengono effettuati da ```performance.py```
    * ```[compound].csv``` sono vari composti che possono essere usati come sink singoli. L'unico che viene usato come source è ```met.csv``` ma nulla vieta di usarne altri
  * ```results``` è la cartella in cui fare produrre i risultati. Il ```.gitignore``` la ignora, per cui non c'è bisogno di preoccuparsi che i risultati vadano in conflitto con quelli degli altri.

## performance.py usage
```
Please specify either --ota, --otc, --oto or --dia
usage: performance.py [-h] [-s] [--min_steps MIN_STEPS] [--max_steps MAX_STEPS]
[--ota] [--otc] [--oto] [--dia] [--min_diam MIN_DIAM] [--max_diam MAX_DIAM] [--out_file OUT_FILE]

Run performance tests for RetroPath2

optional arguments:
  -h, --help            show this help message and exit
  -s, --server          Run tests in server mode
  --min_steps MIN_STEPS
                        Minimum number of steps
  --max_steps MAX_STEPS
                        Maximum number of steps
  --ota                 Run one-to-all test
  --otc                 Run one-to-central test
  --oto                 Run one-to-one test
  --dia                 Run incremental diameters on one-to-one test
  --min_diam MIN_DIAM   Minimum diameter for one-to-one test
  --max_diam MAX_DIAM   Maximum diameter for one-to-one test
  --out_file OUT_FILE   Specify active outputs number
```
Ci sono quattro modalità:
* ota: test "one-to-all" --> usa come source ```met.csv``` e come sink ```byocic_cmpd_rpcleaned.csv```.
* otc: test "one-to-central" --> stessa source (met) ma usa come sink ```gem_central_rp```
* oto: test "one-to-one" --> usa stessa source ma usa un sink diverso (ma sempre un singolo metabolita) a seconda degli step impostati
* dia: test "one-to-one" al variare del diametro delle regole di reazione considerato. Parametri specifici:
  * ```min_diam``` per impostare l'estremo minore del range di diametri da considerare
  * ```max_diam``` estremo massimo
  * ```--min_diam 6 --max_diam 7``` testerà regole con range di diametero 6-6 e poi 6-7

Parametri generali:
* ```s``` per scrivere nei log dei test effettuati sul server. Se non usato scrive nei test effettuati da computer personale
* ```min_steps``` e ```max_steps``` per selezionare il range di steps da testare
  * ```--min_steps 1 --max_steps 4``` testa percorsi di lunghezza 1, poi 2, poi 3 e infine 4


# Usare RetroPath
Github repo del wrapper disponibile [qui](https://github.com/brsynth/RetroPath2-wrapper)

E' necessario attivare l'ambiente virtuale di conda per poter usare retropath. Tutti i file dei composti in ```data_files/``` sono compatibili per essere usati come source/sink.

E' runnabile da command-line: ```python -m retropath2_wrapper <sink_file> <rules_file> <outdir> [-h] [--source_file SOURCE_FILE] [--max_steps MAX_STEPS] [--topx TOPX] [--dmin DMIN] [--dmax DMAX]```
