i due parser creano dei file tsv che possono essere visualizzati aprendo il file cytoscapejs.html, tenendo tutti e tre nella stessa cartella
però bisogna aprire un server virtuale in locale per farlo funzionare e per fargli prendere correttamente i file
per esempio scrivendo: py -m http.server 8000
e andando sul sito: http://localhost:8000/cytoscapejs.html
ancora non funzionano al 100% perchè alcuni legami che sulla mappa di kegg ci sono non li capisce