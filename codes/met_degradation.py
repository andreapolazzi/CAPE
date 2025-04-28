from prova_launch import launch_retropath
from prove_percorsi_noti import get_inchi_from_kegg

# can't obtain pet as is polymer
def obtain_pet():
    kegg_id = "C21451"
    inchi = get_inchi_from_kegg(kegg_id)
    return inchi

def obtain_met():
    return get_inchi_from_kegg("C21450")

def met_to_terephtalate():
    launch_retropath(rules_file="../retrorules_rr02_rp2_flat_forward.csv",
                     source_file="./data_files/met.csv",
                     sink_file="./data_files/terephtalate.csv",
                     pathway_length="1",
                     min_diameter="6",
                     max_diameter="16",
                     max_structures="100",
                     )
    
def met_to_terephtalate_multiple():
    rules_range = (("4", "6"), ("6", "10"), ("10", "14"), ("14", "16"))
    for ranges in rules_range:
        launch_retropath(rules_file="../retrorules_rr02_rp2_flat_forward.csv",
                         source_file="./data_files/met.csv",
                         sink_file="./data_files/terephtalate.csv",
                         pathway_length="1",
                         min_diameter=ranges[0],
                         max_diameter=ranges[1],
                         max_structures="100",
                         output_dir=f"./results/met_ter_range_{ranges[0]}-{ranges[1]}",
                         )
        

def obtain_ter_dihyrodiol():
    return get_inchi_from_kegg("C06318")
    
def met_to_terephtalate_dihydrodiol():
    launch_retropath(rules_file="../retrorules_rr02_rp2_flat_forward.csv",
                     source_file="./data_files/met.csv",
                     sink_file="./data_files/ter_dihydrodiol.csv",
                     pathway_length="2",
                     min_diameter="6",
                     max_diameter="16",
                     max_structures="100",
                     )
    
def obtain_dihydroxy_benzoate():
    return get_inchi_from_kegg("C00230")

def met_to_dihydroxy_benzoate():
    launch_retropath(rules_file="../retrorules_rr02_rp2_flat_forward.csv",
                     source_file="./data_files/met.csv",
                     sink_file="./data_files/dihydroxy_benzoate.csv",
                     pathway_length="3",
                     min_diameter="6",
                     max_diameter="16",
                     max_structures="100",
                     )


def obtain_beta_muconate():
    return get_inchi_from_kegg("C01163")

def met_to_beta_muconate():
    launch_retropath(rules_file="../retrorules_rr02_rp2_flat_forward.csv",
                     source_file="./data_files/met.csv",
                     sink_file="./data_files/beta_carboxy_muconate.csv",
                     pathway_length="4",
                     min_diameter="6",
                     max_diameter="16",
                     max_structures="100",
                     )
    

def obtain_muconolactone():
    return get_inchi_from_kegg("C01278")

def met_to_muconolactone():
    launch_retropath(rules_file="../retrorules_rr02_rp2_flat_forward.csv",
                     source_file="./data_files/met.csv",
                     sink_file="./data_files/muconolactone.csv",
                     pathway_length="5",
                     min_diameter="6",
                     max_diameter="16",
                     max_structures="100",
                     #max_structures="50",
                     )
    




if __name__ == "__main__":
    #print(obtain_met())
    #met_to_terephtalate()
    print("MET to Terephtalate correct")
    #met_to_terephtalate_multiple()
    print("MET to Terephtalate multiple identified rules for diameters 6, 14, 16 (check) all with EC:3.1.1.102")
    print("The results (1 everytime) are not the same as running all diameters together\n")
    
    #print(obtain_ter_dihyrodiol())
    #met_to_terephtalate_dihydrodiol()
    print("MET to Terephthalate-1,2-cis-dihydrodiol correct")
    print("Produced 1 result\n")

    #print(obtain_dihydroxy_benzoate())
    #met_to_dihydroxy_benzoate()
    print("MET to 3,4-dihydroxybenzoate correct")
    print("Only checked presence of ec\n")

    #print(obtain_beta_muconate())
    #met_to_beta_muconate()
    print("MET to beta-Carboxy-cis,cis-muconate correct")
    print("Only checked presence of ec\n")

    #print(obtain_muconolactone())
    met_to_muconolactone()
    print("Couldn't complete 5th step w/ 100 structures, probably due to memory or time")
    print("W/ 50 an error occured")