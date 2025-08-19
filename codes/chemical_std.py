from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import inchi

def standardize(inchi_str, 
                standardize=True,
                neutralize=True, 
                strip_stereo=False, 
                canonicalize_tautomers=False):
    """
    Standardize molecule to KEGG-like representation.
    """
    # Load molecule from InChI
    mol = Chem.MolFromInchi(inchi_str)
    if mol is None:
        return None
    
    #0. Standardize to avoid invalid groups
    if standardize:
        mol = rdMolStandardize.Normalize(mol)
    
    # 1. Neutralize charges if requested
    if neutralize:
        mol = rdMolStandardize.Uncharger().uncharge(mol)  
        #mol = uncharger.uncharge(mol)
    
    # 2. Canonicalize tautomers
    if canonicalize_tautomers:
        #tautomer_enumerator = rdMolStandardize.TautomerEnumerator()
        #mol = tautomer_enumerator.Canonicalize(mol)
        mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
    
    # 3. Strip stereochemistry if requested
    if strip_stereo:
        Chem.RemoveStereochemistry(mol)
    
    # 4. Get standardized InChI and InChIKey
    std_inchi = inchi.MolToInchi(mol)
    std_inchikey = inchi.InchiToInchiKey(std_inchi)
    
    return {"inchi": std_inchi, "inchikey": std_inchikey, "mol": mol}


if __name__ == '__main__':
    # Example: L-lactic acid InChI (charged vs neutral)
    lactic_inchi = "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)"
    lactate_inchi = "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/p-1"
    Dlactic_inchi = "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m1/s1"
    Llactic_inchi = "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m0/s1"
    pyruvate_inchi = "InChI=1S/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)/p-1"

    # standardize it in KEGG-like style
    standardized_lactic = standardize(lactic_inchi, 
                                    neutralize=True, 
                                    strip_stereo=False, 
                                    canonicalize_tautomers=False)
    standardized_lactate = standardize(lactate_inchi, 
                                    neutralize=True, 
                                    strip_stereo=False, 
                                    canonicalize_tautomers=False)
    standardized_Dlactic = standardize(Dlactic_inchi, 
                                    neutralize=True, 
                                    strip_stereo=False, 
                                    canonicalize_tautomers=False)

    print("standardized lactic InChI:", standardized_lactic["inchi"])
    print("standardized lactic InChIKey:", standardized_lactic["inchikey"])
    print("standardized lactate InChI:", standardized_lactate["inchi"])
    print("standardized lactate InChIKey:", standardized_lactate["inchikey"])
    print("standardized pyruvate InChI:", standardized_Dlactic["inchi"])
    print("standardized pyruvate InChIKey:", standardized_Dlactic["inchikey"])

    m = Chem.MolFromInchi(Llactic_inchi)
    m2 = rdMolStandardize.Normalize(m)
    m2_inchi = inchi.MolToInchi(m2)
    m2_inchikey = inchi.InchiToInchiKey(m2_inchi)
    print(m2_inchi, m2_inchikey)
