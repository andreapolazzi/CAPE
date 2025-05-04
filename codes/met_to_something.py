from prova_launch import launch_retropath

def one_step():
    launch_retropath(rules_file="../retrorules_rr02_rp2_flat_forward.csv",
                     source_file="./data_files/met.csv",
                     sink_file="./data_files/biocyc_cmpd_rpcleaned.csv",
                     pathway_length="1",
                     min_diameter="6",
                     max_diameter="16",
                     max_structures="100",
                     )
    
if __name__ == "__main__":
    
    one_step()