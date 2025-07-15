from prova_launch import launch_retropath
from results_parser import find_pathways, results_to_df

def n_step_met(n: int):
    launch_retropath(rules_file="../retrorules_rr02_rp2_flat_forward.csv",
                     source_file="./data_files/met.csv",
                     sink_file="./data_files/gem_central2_rp.csv",
                     pathway_length=f"{n}",
                     min_diameter="6",
                     max_diameter="16",
                     max_structures="100",
                     )

    
if __name__ == "__main__":
    
    #n_step_met(1)
    print("One-step long pathway produced results. None in sink. Reasonable products")

    #n_step_met(2)
    print("Two-step long pathway produced results. Many in sink. Check results better. " \
    "At least 1 pathway is consistent with annotated pathway")

    #n_step_met(3)
    print("Three-step long pathway produced results. Many in sink. Check results better. " \
    "At least 1 pathway is consistent with annotated pathway")

    #n_step_met(4)
    print("One-step long pathway produced results. Strange output: very few results. Muconate not found")

    # df = results_to_df("results/results.csv")
    # print(df)
    # paths = find_pathways(df, '3-4-DIHYDROXYBENZOATE')
    # print(paths)

    #n_step_met(4)
    print("Nothing interesting")

    n_step_met(3)