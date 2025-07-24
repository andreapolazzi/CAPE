from rp_classes import RP

return_code, stdout, stderr = RP.rp_single_source(
    sink="data_files/kegg_rp.csv", 
    name="O-Xylene",
    inchi="InChI=1S/C8H10/c1-7-5-3-4-6-8(7)2/h3-6H,1-2H3",
    rules="../retrorules_rr02_rp2_flat_forward.csv",
    outdir="results/",
    timeout=1,
    steps=1
    )
print(return_code, '\n\n', stdout, '\n\n', stderr)