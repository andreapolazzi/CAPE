def primo(s1,s2):
    a=[s1,s2,s1+s2]
    return a, ['ec1','ec2','ec3']


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

