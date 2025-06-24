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
  * ```data_files/``` contiene i file di dati, come i vari sink e i source
  * ```results``` è la cartella in cui fare produrre i risultati. Il ```.gitignore``` la ignora, per cui non c'è bisogno di preoccuparsi che i risultati vadano in conflitto con quelli degli altri.