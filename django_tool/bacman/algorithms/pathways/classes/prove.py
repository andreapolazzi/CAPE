from chem_utils import *
from rp_classes import *


print("=== RETROPATH PIPELINE EXAMPLE ===")

print("Generating the tsv file from RetroPath results...")
RpOutput.tsv_generator(results_path="/work/igem/iGEM/django_tool/bacman/algorithms/pathways/rp_results/results.csv",
					   outdir="/work/igem/iGEM/django_tool/bacman/algorithms/pathways/rp_results/")

print("Standardizing the edges...")
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