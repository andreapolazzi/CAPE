from chem_utils import *
from rp_classes import *


print("=== RETROPATH PIPELINE EXAMPLE ===")

print("Generating retropath results...")
RP.django_call(name="Methyl_tert_butyl_ether", 
			   inchi="InChI=1S/C5H12O/c1-5(2,3)6-4/h1-4H3")

print("Generating the tsv file from RetroPath results...")
RpOutput.tsv_generator(results_path="/work/igem/iGEM/django_tool/bacman/algorithms/pathways/rp_results/results.csv",
					   outdir="/work/igem/iGEM/django_tool/bacman/algorithms/pathways/rp_results/")

print("Standardizing the edges...")
# EXPECTS A SINGLE UNIQUE KEGG FILE FOR THE EDGES
Std.main(tsv_edges_path="/work/igem/iGEM/django_tool/bacman/algorithms/pathways/rp_results/rp_edges.tsv",
		 edges_std_path="/work/igem/iGEM/django_tool/bacman/algorithms/pathways/rp_results/rp_std_edges.tsv",
		 processes=64)

print("Adding weights to the edges...")
Penalty.main(std_path="/work/igem/iGEM/django_tool/bacman/algorithms/pathways/rp_results/rp_std_edges.tsv",
			 weighted_path="/work/igem/iGEM/django_tool/bacman/algorithms/pathways/rp_results/rp_weighted_edges.tsv")

# -----------------------------------

print("=== KEGG PIPELINE EXAMPLE ===")

print("Standardizing the edges...")
Std.main(tsv_edges_path="/work/igem/standardization/roa_inchi_edges.tsv",
		 edges_std_path="/work/igem/graphs/kegg_graph/roa_std_edges2.tsv",
		 processes=64)

print("Adding weights to the edges...")
Penalty.main(std_path="/work/igem/graphs/kegg_graph/roa_std_edges2.tsv",
			 weighted_path="/work/igem/graphs/kegg_graph/roa_weighted_edges2.tsv")

# -----------------------------------

print("=== UNIQUE PIPELINE FOR BOTH ===")

# first make RP output compatible with KEGG data
print("Generating the tsv file from RetroPath results...")
RpOutput.tsv_generator(results_path="/work/igem/iGEM/django_tool/bacman/algorithms/pathways/rp_results/results.csv",
					   outdir="/work/igem/iGEM/django_tool/bacman/algorithms/pathways/rp_results/")

# COMBINE KEGG AND RP DATA
print("Combining KEGG and RetroPath data...")
# function to do it

print("Standardizing the edges...")
Std.main(tsv_edges_path="/path/to/combined/edges.tsv",  # hypothetical
		 edges_std_path="/path/to/combined/std_edges.tsv",  # hypothetical
		 processes=64)

print("Adding weights to the edges...")
Penalty.main(std_path="/path/to/combined/std_edges.tsv",  # hypothetical
			 weighted_path="/path/to/combined/weighted_edges.tsv")  # hypothetical