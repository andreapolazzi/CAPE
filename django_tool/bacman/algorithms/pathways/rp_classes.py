import subprocess
import pandas as pd
from pubchempy import get_compounds
import requests

class RP():

    
    CONDA_INIT = "source /work/igem/anaconda3/etc/profile.d/conda.sh"
    ENV_NAME = "retropath_env"


    @staticmethod
    def rp_source_file_noenv(sink: str, rules:str, outdir:str,
                       source: str,
                       steps=3, structures=100, dmin=6, dmax=16, timeout=60):
        """ Launches retropath2 using a source file configuration """
        command = [
            "python", "-m", "retropath2_wrapper",
            sink,
            rules,
            outdir,
            "--source_file", source,
            "--max_steps", str(steps),
            "--topx", str(structures),
            "--dmin", str(dmin),
            "--dmax", str(dmax),
            "--msc_timeout", str(timeout)  # minutes
        ]

        try:
            # Print command for debugging
            # print("Executing command:", ' '.join(command))
            
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=False  # Don't raise exception on non-zero return code
            )
            return result.returncode, result.stdout, result.stderr
            #return result.returncode
            
        except Exception as e:
            return 99, "", str(e)
        
    @staticmethod
    def rp_single_source_noenv(sink: str, rules:str, outdir:str, 
                         name: str, inchi: str,
                       steps=3, structures=100, dmin=6, dmax=16, timeout=60):
        """ Launches retropath2 using source name and InChI """
        command = [
            "python", "-m", "retropath2_wrapper",
            sink,
            rules,
            outdir,
            "--source_name", name,
            "--source_inchi", inchi,
            "--max_steps", str(steps),
            "--topx", str(structures),
            "--dmin", str(dmin),
            "--dmax", str(dmax),
            "--msc_timeout", str(timeout)  # minutes
        ]
    
        try:
            # Print command for debugging
            # print("Executing command:", ' '.join(command))
            
            # Run the command using subprocess.run
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=False  # Don't raise exception on non-zero return code
            )
            
            return result.returncode, result.stdout, result.stderr
            
        except Exception as e:
            return 99, "", str(e)
        

    @staticmethod
    def rp_source_file(sink: str, rules: str, outdir: str,
                    source: str,
                    steps=3, structures=100, dmin=6, dmax=16, timeout=60):
        """ 
        Launches retropath2 using a source file configuration inside retropath_env 

        Return codes:
            0: No error
            1: File is not found
            2: Running the RetroPath2.0 Knime program produced an OSError
            3: InChI is malformated
            10: Source has been found in the sink (warning)
            11: No solution is found (warning)
            99: unexpected error (not defined in retropath module)
        """
        command = f"""
        {RP.CONDA_INIT} && conda activate {RP.ENV_NAME} && \
        python -m retropath2_wrapper {sink} {rules} {outdir} \
        --source_file {source} \
        --max_steps {steps} \
        --topx {structures} \
        --dmin {dmin} \
        --dmax {dmax} \
        --msc_timeout {timeout}
        """
        try:
            result = subprocess.run(command, shell=True, executable="/bin/bash",
                                    capture_output=True, text=True, check=False)
            return result.returncode, result.stdout, result.stderr
        except Exception as e:
            return 99, "", str(e)

    @staticmethod
    def rp_single_source(sink: str, rules: str, outdir: str,
                        name: str, inchi: str,
                        steps=3, structures=100, dmin=6, dmax=16, timeout=60):
        """ 
        Launches retropath2 using source name and InChI inside retropath_env 
        
        Return codes:
            0: No error
            1: File is not found
            2: Running the RetroPath2.0 Knime program produced an OSError
            3: InChI is malformated
            10: Source has been found in the sink (warning)
            11: No solution is found (warning)
            99: unexpected error (not defined in retropath module)
        """
        command = f"""
        {RP.CONDA_INIT} && conda activate {RP.ENV_NAME} && \
        python -m retropath2_wrapper {sink} {rules} {outdir} \
        --source_name "{name}" \
        --source_inchi "{inchi}" \
        --max_steps {steps} \
        --topx {structures} \
        --dmin {dmin} \
        --dmax {dmax} \
        --msc_timeout {timeout}
        """
        try:
            result = subprocess.run(command, shell=True, executable="/bin/bash",
                                    capture_output=True, text=True, check=False)
            return result.returncode, result.stdout, result.stderr
        except Exception as e:
            return 99, "", str(e)


class RpGraph():

    @staticmethod
    def load_res(file_path):
        """ Load the results from a CSV file into a pandas DataFrame. """
        return pd.read_csv(file_path)
    
    @staticmethod
    def get_name(inchi):
        """ Get the name of a compound from its InChI """
        if inchi == 'InChI=1S':
            return 'H+'
        compounds = get_compounds(inchi, 'inchi')
        if compounds:
            if compounds[0].synonyms:
                return compounds[0].synonyms[0]
            elif compounds[0].iupac_name:
                return compounds[0].iupac_name
            else:
                return 'Unknown'
        return 'Unknown'
    
    @staticmethod
    def add_names(df):
        """
        Add 'Substrate Name' and 'Product Name' columns to the DataFrame based on InChI.
        Queries pubchem to get names.
        """

        df['Substrate Name'] = df['Substrate InChI'].apply(RpGraph.get_name)
        df['Product Name'] = df['Product InChI'].apply(RpGraph.get_name)
        
        return df
    
    @staticmethod
    def nodes_to_tsv(df, outdir):
        """ Save the unique substrate and product names to a TSV file """
        
        unique_substrates = pd.Series(df['Substrate InChI'].unique())
        unique_products = pd.Series(df['Product InChI'].unique())

        all_compounds = pd.concat([unique_substrates, unique_products]).unique()
        all_compounds = pd.DataFrame({'inchi': all_compounds})
        all_compounds['name'] = all_compounds['inchi'].apply(RpGraph.get_name)
        all_compounds = all_compounds[['name', 'inchi']]

        all_compounds.to_csv(outdir+'rp_nodes.tsv', sep='\t', index=False, header=True)

    @staticmethod
    def ec_to_name(ec: str):
        """ Convert EC number (format "x.x.x.x", not "EC:x.x.x.x") to a human-readable name. """
        try:
            response = requests.get(f'https://rest.kegg.jp/get/ec:{ec}')
            if response.status_code == 200:
                lines = response.text.split('\n')
                for line in lines:
                    if line.startswith('NAME'):
                        return line.split('NAME')[1].strip().strip(';')
            return 'Unknown'
        except:
            return 'Unknown'
        
    @staticmethod
    def edges_to_tsv(df, outdir):
        """ Save the edges (reactions) to a TSV file. """
        df['EC number'] = df['EC number'].apply(lambda x: x.strip('[]'))
        #print(df['EC number'])
        edges_df = df[['EC number', 'Substrate InChI', 'Product InChI']]
        #print(edges_df)
        edges_df['enzyme_name'] = edges_df['EC number'].apply(lambda x: x if x == 'NOEC' else RpGraph.ec_to_name(x))
        edges_df['color'] = 'yellow'
        edges_df = edges_df.rename(columns={'EC number': 'EC', 'Substrate InChI': 'source', 'Product InChI': 'target'})
        edges_df = edges_df[['enzyme_name', 'EC', 'source', 'target', 'color']]
        edges_df.to_csv(outdir+'rp_edges.tsv', sep='\t', index=False, header=True)
        #print(edges_df)


    @staticmethod
    def tsv_generator(file_path, outdir):
        """ Generates nodes and edges tsv (2 separate files) """
        
        df = RpGraph.load_res(file_path)      

        RpGraph.nodes_to_tsv(df, outdir)
        print("Nodes saved to TSV file.")

        RpGraph.edges_to_tsv(df, outdir)
        print("Edges saved to TSV file.")