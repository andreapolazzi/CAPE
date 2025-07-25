def primo(s1,s2):
    a=[ 'Path #1 (total length: 7)\n 1. Mono(ethylene terephtalate) ---[3.1.1.102]---> Terephtalate (KEGG white edge)\n 2. Terephtalate ---[1.14.12.15]---> Terephtalate-1,2-cis-dihydrodiol (KEGG white edge)\n 3. Terephtalate-1,2-cis-dihydrodiol ---[1.3.1.53]---> 3,4-Dihydroxy-benzoate (KEGG white edge)\n 4. 3,4-Dihydroxy-benzoate ---[1.13.11.3]---> beta-Carboxy-muconate (KEGG green edge)\n 5. beta-Carboxy-muconate ---[5.5.1.2]---> gamma-Carboxy-muconolactone (KEGG green edge)\n 6. gamma-Carboxy-muconolactone ---[4.1.1.44]---> 3-Oxoadipate-enol-lactone (KEGG green edge)\n 7. 3-Oxoadipate-enol-lactone ---[3.1.1.24]---> 3-Oxo-adipate (KEGG green edge)',
        'Path #2 (total length: )\n 1. Mono(ethylene terephtalate) ---[3.1.1.102]---> Terephtalate (KEGG white edge)\n 2. Terephtalate ---[1.14.12.15]---> Terephtalate-1,2-cis-dihydrodiol (KEGG white edge)\n 3. Terephtalate-1,2-cis-dihydrodiol ---[1.3.1.53]---> 3,4-Dihydroxy-benzoate (KEGG white edge)\n 4. 3,4-Dihydroxy-benzoate ---[1.14.13.64]---> Hydroxyquinol (KEGG white edge)\n 5. Hydroxyquinol ---[1.13.11.37]---> 3-Hydroxy-cis-cis-muconate (KEGG green edge)\n 6. 3-Hydroxy-cis-cis-muconate ---> 2-Maleylacetate \n 7. 2-Maleylacetate ---[1.3.1.32]---> 3-Oxo-adipate (KEGG green edge)',
        'No path 3'
      ]

    b=[['3.1.1.102', '1.14.12.15', '1.3.1.53'], ['3.1.1.102', '1.14.12.15', '1.3.1.53','1.14.13.64'],[]]
    return a,b


def is_curated(pollutant) -> bool:
    """ Returns True if the pollutant is present in KEGG paths """
    pass

def path_finder(pollutant, target, graph):
    """ Returns the path in the graph from pollutant to target """
    pass

def retropath_workflow(pollutant, poll_inchi):
    """ Runs RetroPath and produces tsv for nodes and edges """
    pass

def merge_graphs():
    """ Creates one single graph containing info from RetroPath and KEGG """
    pass


def general_workflow(pollutant: str, poll_inchi: str, target: str):
    """ General workflow of Step One: Paths Finding"""

    if is_curated(pollutant):
        path_finder()

    else:
        retropath_workflow(pollutant, poll_inchi)
        merge_graphs()
        path_finder()

