from rp_classes import RP

return_code, stdout, stderr = RP.rp_single_source(
    name="bp3",
    inchi="InChI=1S/C14H12O3/c1-17-11-7-8-12(13(15)9-11)14(16)10-5-3-2-4-6-10/h2-9,15H,1H3",
    sink="data_files/kegg_rp.csv",
    rules="../retrorules_rr02_rp2_flat_forward.csv",
    outdir="results/",
    timeout=1,
    steps=1
    )
print(return_code, '\n\n', stdout, '\n\n', stderr)