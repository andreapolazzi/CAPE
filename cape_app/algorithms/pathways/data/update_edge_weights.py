import os
import csv
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs,rdFingerprintGenerator
import numpy as np


def transf_distance(substrate_id: str, product_id: str, format="inchi") -> float:
    """Tanimoto-based chemical distance (1-similarity)."""
    if not substrate_id or not product_id:
        return None
    mol1 = Chem.MolFromSmiles(substrate_id) if format=="smiles" else Chem.MolFromInchi(substrate_id)
    mol2 = Chem.MolFromSmiles(product_id) if format=="smiles" else Chem.MolFromInchi(product_id)
    if mol1 is None or mol2 is None:
        return None
    try:
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        fp1 = gen.GetFingerprint(mol1)
        fp2 = gen.GetFingerprint(mol2)
        return 1 - DataStructs.TanimotoSimilarity(fp1, fp2)
    except:
        return None


def dist_penalty(distance: float) -> float:
    """Normalize distance into a 0-5 penalty."""
    if distance is None:
        return 0
    if distance <= 0.70:
        penalty = distance * 2.0
    elif distance < 0.80:
        penalty = 0.60 * 2.0 + (distance - 0.60) * 5.0
    else:
        penalty = 0.60 * 2.0 + 0.20 * 5.0 + (distance - 0.80) * 12.0
    return penalty * 5 / 4.6


def mw_change(substrate_mw: float, product_mw: float) -> float:
    """Molecular weight change penalty (0-25)."""
    if substrate_mw is None or product_mw is None or substrate_mw <= 0:
        return 0
    mw_diff = substrate_mw - product_mw
    percent_change = abs(mw_diff) / substrate_mw * 100
    if mw_diff >= 0:
        penalty = 3.0 * (1 - np.exp(-percent_change / 40))
    else:
        penalty = 3.0 * (1 - np.exp(-percent_change / 25))
    return (min(penalty, 3.0)*5/3)**2


def transf_penalty(rp_score: float, distance_penalty: float, mw_penalty: float) -> float:
    """Combine penalties into final edge weight."""
    if rp_score is None:
        rp_score = 0.0
    if distance_penalty is None:
        distance_penalty = 0.0
    if mw_penalty is None:
        mw_penalty = 0.0
    return rp_score + mw_penalty #+ distance_penalty (excluded for now)


blacklist = [
    # Energy carriers / Electron carriers
    "ZKHQWZAMYRWXGA-KQYNXXCUSA-N", # ATP
    "XTWYTFMLZFPYCI-KQYNXXCUSA-N", # ADP
    "UDMBCSSLTHHNCD-KQYNXXCUSA-N", # AMP
    "XKMLYUALXHKNFT-UUOKFMHZSA-N", # GTP
    "QGWNDRXFNXRZMB-UUOKFMHZSA-N", # GDP
    "RQFCJASXJCIDSX-UUOKFMHZSA-N", # GMP
    "PGAVKCOVUIYSFO-XVFCMESISA-N", # UTP
    "XCCTYIAWTASOJW-XVFCMESISA-N", # UDP
    "DJJCXFVJDGTHFX-XVFCMESISA-N", # UMP
    "PCDQPRRSZKQHHS-XVFCMESISA-N", # CTP
    "ZWIADYZPOWUWEW-XVFCMESISA-N", # CDP
    "IERHLVCPSMICTF-XVFCMESISA-N", # CMP
    "NBIIXXVUZAFLBC-UHFFFAOYSA-N", # Orthophosphate
    "XPPKVPWEQAFLFU-UHFFFAOYSA-N", # Diphosphate
    "UNXRWKVEANCORM-UHFFFAOYSA-N", # Polyphosphate
    "YTNIXZGTHTVJBW-SCRDCRAPSA-N", # Reduced FMN
    "BAWFJGJZGIEFAR-NNYOXOHSSA-O", # NAD+
    "BOPGDPNILDQYTO-NNYOXOHSSA-N", # NADH
    "XJLXINKUBYWONI-NNYOXOHSSA-O", # NADP+
    "ACFIXJIJDZMPPO-NNYOXOHSSA-N", # NADPH
    "VWWQXMAJTJZDQX-UYBVJOGSSA-N", # FAD
    "YPZRHBJKEMOYQH-UYBVJOGSSA-N", # FADH2
    "FVTCRASFADXXNN-SCRDCRAPSA-N", # FMN
    # Coenzymes and group carriers
    "RGJOEKWQDUBAIZ-DRCCLKDXSA-N", # CoA, Coenzyme A
    "ZSLZBFCDCINBPY-KMYLAXNMSA-N", # Acetyl-CoA
    "MSTNYGQPCMXVAQ-PXYINDEMSA-N", # Tetrahydrofolate
    "YBJHBAHKTGYVGT-ZKWXMUAHSA-N", # Biotin
    "NGVDGCNFYWLIFO-UHFFFAOYSA-N", # Pyridoxal phosphate
    "AYEKOFBPNLCAJY-UHFFFAOYSA-O", # Thiamine Diphosphate
    # Gases and small byproducts
    "CURLTUGMZLYLDI-UHFFFAOYSA-N", # CO2, Carbon Dioxide
    "XLYOFNOQVPJJNP-UHFFFAOYSA-N", # H2O, Water
    "QGZKDVFQNNGYKY-UHFFFAOYSA-N", # NH3, Ammonia
    "GPRLSGONYQIRFK-UHFFFAOYSA-N", # H+, Hydrogen
    "MYMOFIZGZYHOMD-UHFFFAOYSA-N", # O2, Oxygen
    "MHAJPDPJQMAIIY-UHFFFAOYSA-N", # H2O2, Hydrogen Peroxide
    "VNWKTOKETHGBQD-UHFFFAOYSA-N", # Methane, CH4
    "BDAGIHXWWSANSR-UHFFFAOYSA-N", # Formate
    "WSFSSNUMVMOOMR-UHFFFAOYSA-N", # Formaldehyde
    # Sugars and nucleotides in transfer reactions
    "HSCJRCZFDFQWRP-JZMIEXBBSA-N", # UDP-Glucose
    "RZZPDXZPRHQOCG-OJAKKHQRSA-N", # CDP-Choline
    "JSRLJPSBLDHEIO-SHYZEUOFSA-N", # dUMP
    "GYOZYWVXFNDGLU-XLPZGREQSA-N", # dTMP
    # Miscellaneous
    "RWSXRVCMGQZWBV-UHFFFAOYSA-N", # Glutathione
    "PHIQHXFUZVPYII-UHFFFAOYSA-O", # Carnitine
    "CVSVTCORWBXHQV-UHFFFAOYSA-N", # Creatine
    "DRBBFCLWYRJSJZ-UHFFFAOYSA-N", # Phosphocreatine
    "CRADQHJNYKKORV-UHFFFAOYSA-N", # Ubiquinone
    "RJUKAZDREBHPDA-UHFFFAOYSA-N", # Ubiquinol
    "GAOZTHIDHYLHMS-UHFFFAOYSA-N", # Ferricytochrome b5
    "LJZLHDYZFURWEE-UHFFFAOYSA-N", # Ferricytochrome cL
]

def add_weights(pathway_folder: str):
    dist_cache = {}
    for file_name in os.listdir(pathway_folder):
        if not file_name.endswith("_edges.tsv"):
            continue
        file_path = os.path.join(pathway_folder, file_name)
        temp_file = file_path + ".tmp"
        
        with open(file_path, newline="", encoding="utf-8") as f_in, \
             open(temp_file, "w", newline="", encoding="utf-8") as f_out:
            
            reader = csv.DictReader(f_in, delimiter="\t")
            fieldnames = reader.fieldnames.copy()
            if "weight" not in fieldnames:
                fieldnames.append("weight")
            writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            
            for row in reader:
                # Blacklist check
                if row.get("source_inchi_key") in blacklist or \
                row.get("target_inchi_key") in blacklist or \
                row.get("source_mol_weight") == '999999.0' or \
                row.get("target_mol_weight") == '999999.0':
                    row["weight"] = 999999
                    writer.writerow(row)
                    continue
                
                # Chemical distance
                source_inchi = row.get("source_inchi")
                target_inchi = row.get("target_inchi")
                key = tuple(sorted((source_inchi, target_inchi)))
                if key in dist_cache:
                    dist = dist_cache[key]
                else:
                    dist = transf_distance(source_inchi, target_inchi)
                    dist_cache[key] = dist
                d_penalty = dist_penalty(dist)
                
                # MW penalty
                try:
                    s_mw = float(row.get("source_mol_weight",0))
                except:
                    s_mw = 0.0
                try:
                    t_mw = float(row.get("target_mol_weight",0))
                except:
                    t_mw = 0.0
                mw_pen = mw_change(s_mw, t_mw)
                
                # RP score
                try:
                    rp_score = float(row.get("RP_score",0))
                except:
                    rp_score = 0.0
                
                # Final edge weight
                col = row.get("color")
                penalty = transf_penalty(rp_score, d_penalty, mw_pen)
                if col == "green":
                    row["weight"] = penalty
                elif col == "white":
                    row["weight"] = penalty*2
                else:
                    row["weight"] = penalty*3
                writer.writerow(row)

        # Replace original file
        os.replace(temp_file, file_path)


if __name__ == "__main__":

    add_weights("./Pathways")

