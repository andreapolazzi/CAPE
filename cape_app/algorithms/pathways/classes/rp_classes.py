import subprocess
import pandas as pd
from pubchempy import get_compounds
import requests
import os
import shutil
import os
from pathlib import Path

class RP():

    
    CONDA_INIT = "source /work/igem/anaconda3/etc/profile.d/conda.sh"  # Adjust this path to your conda installation
    ENV_NAME = "retropath_env"  # Name of the conda environment with retropath installed


    @staticmethod
    def rp_source_file_noenv(sink: str, rules:str, outdir:str,
                       source: str,
                       steps=3, structures=100, dmin=6, dmax=16, timeout=60) -> tuple:
        """
        Launches retropath2 using source file
        Requires a retropath conda environment activated prior its call

        Args:
            For more info on RetroPath parameters, see their documentation
            sink: path to sink file
            rules: path to rules file
            outdir: path to the folder where to save results
            source: path to the source file. Must contain:
                name of the source compound (no need to be precise)
                inchi: inchi code of the source compound (needs to be exact). Must start with "InChI="
            steps: max iterations RetroPath is allowed to compute
            structures: max number of structures kept after each iteration
            dmin: min diameter for the reaction rules
            dmax: max diameter for the reaction rules
            timeout: max time RetroPath is allowed to run, in minutes
        
        Returns:
            result.returncode: return code
            result.stdout: captured retropath text output
            result.stderr: exception raised or failed RetroPath call
        
        Return codes:
            0: No error
            1: File is not found
            2: Running the RetroPath2.0 Knime program produced an OSError
            3: InChI is malformated
            10: Source has been found in the sink (warning)
            11: No solution is found (warning)
            99: unexpected error (not defined in retropath module)
        
        """
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
    def rp_single_source_noenv(sink: str, rules:str, outdir:str,
                               name: str, inchi: str,
                               steps=3, structures=100, dmin=6, dmax=16, timeout=60) -> tuple:
        """
        Launches retropath2 using source name and InChI
        Requires a retropath conda environment to be activated prior its call

        Args:
            For more info on RetroPath parameters, see their documentation
            sink: path to sink file
            rules: path to rules file
            outdir: path to the folder where to save results
            name: name of the source compound (no need to be precise)
            inchi: inchi code of the source compound (needs to be exact)
                    Must start with "InChI="
            steps: max iterations RetroPath is allowed to compute
            structures: max number of structures kept after each iteration
            dmin: min diameter for the reaction rules
            dmax: max diameter for the reaction rules
            timeout: max time RetroPath is allowed to run, in minutes
        
        Returns:
            result.returncode: return code
            result.stdout: captured retropath text output
            result.stderr: exception raised or failed RetroPath call
        
        Return codes:
            0: No error
            1: File is not found
            2: Running the RetroPath2.0 Knime program produced an OSError
            3: InChI is malformated
            10: Source has been found in the sink (warning)
            11: No solution is found (warning)
            99: unexpected error (not defined in retropath module)
        
        """
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
            # Run the command
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
                    steps=3, structures=100, dmin=6, dmax=16, timeout=60) -> tuple:
        """ 
        Launches retropath2 using source file
        Requires a retropath conda environment available to be activated

        Args:
            For more info on RetroPath parameters, see their documentation
            sink: path to sink file
            rules: path to rules file
            outdir: path to the folder where to save results
            source: path to the source file. Must contain:
                name of the source compound (no need to be precise)
                inchi: inchi code of the source compound (needs to be exact). Must start with "InChI="
            steps: max iterations RetroPath is allowed to compute
            structures: max number of structures kept after each iteration
            dmin: min diameter for the reaction rules
            dmax: max diameter for the reaction rules
            timeout: max time RetroPath is allowed to run, in minutes
        
        Returns:
            result.returncode: return code
            result.stdout: captured retropath text output
            result.stderr: exception raised or failed RetroPath call
        
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
                        steps=5, structures=100, dmin=6, dmax=16, timeout=60) -> tuple:
        """ 
        Launches retropath2 using source name and InChI
        Requires a retropath conda environment available to be activated

        Args:
            For more info on RetroPath parameters, see their documentation
            sink: path to sink file
            rules: path to rules file
            outdir: path to the folder where to save results
            name: name of the source compound (no need to be precise)
            inchi: inchi code of the source compound (needs to be exact)
                    Must start with "InChI="
            steps: max iterations RetroPath is allowed to compute
            structures: max number of structures kept after each iteration
            dmin: min diameter for the reaction rules
            dmax: max diameter for the reaction rules
            timeout: max time RetroPath is allowed to run, in minutes
        
        Returns:
            result.returncode: return code
            result.stdout: captured retropath text output
            result.stderr: exception raised or failed RetroPath call
        
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
        
    def clean_results(outdir: str):
        """ Cleans the output directory before running retropath """
        if os.path.exists(outdir):
            shutil.rmtree(outdir)
        os.makedirs(outdir, exist_ok=True)

    def django_call(name: str, inchi: str,
                    sink="./cape_app/algorithms/pathways/data/gsm_gapfilled_sink.csv",
                    rules="./cape_app/algorithms/pathways/data/retrorules_rr02_rp2_flat_forward.csv", 
                    #outdir="./cape_app/algorithms/pathways/rp_results/",
                    outdir="./cape_app/pipeline_results/rp_results",
                    steps=5, structures=100, dmin=6, dmax=16, timeout=60) -> tuple:
        """ 
        Utility to call retropath from django
        Acts a main for retropath execution with single source

        !! WARNING: this will delete all the content of outdir !!

        Args:
            name: name of the source compound (no need to be precise)
            inchi: inchi code of the source compound (needs to be exact)
                    Must start with "InChI="
            sink: path to sink file
            rules: path to rules file
            outdir: path to the folder where to save results
            steps: max iterations RetroPath is allowed to compute
            structures: max number of structures kept after each iteration
            dmin: min diameter for the reaction rules
            dmax: max diameter for the reaction rules
            timeout: max time RetroPath is allowed to run, in minutes
        
        Returns:
            code: return codes from retropath execution
        """
        RP.clean_results(outdir)
        codes = RP.rp_single_source_noenv(sink=sink, rules=rules, outdir=outdir,
                                   name=name, inchi=inchi,
                                   steps=steps, structures=structures,
                                   dmin=dmin, dmax=dmax, timeout=timeout)
        return codes


class RpOutput():

    @staticmethod
    def load_res(file_path) -> pd.DataFrame:
        """ 
        Load the results from a CSV file into a pandas DataFrame.
        
        Args:
            file_path (str): Path to the CSV file.
        Returns:
            pd.DataFrame: DataFrame containing the CSV data.
        """
        return pd.read_csv(file_path)
    
    @staticmethod
    def get_name(inchi: str) -> str:
        """ 
        Get the name of a compound from its InChI
        Queries pubchem to get names.

        Args:
            inchi (str): InChI string of the compound.
        Returns:
            str: Name of the compound or 'Unknown' if not found.
        """
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
    def add_names(df: pd.DataFrame) -> pd.DataFrame:
        """
        Add 'Substrate Name' and 'Product Name' columns to the DataFrame based on InChI.
        Queries pubchem to get names.

        Args:
            df (pd.DataFrame): DataFrame containing 'Substrate InChI' and 'Product InChI' columns.
        Returns:
            pd.DataFrame: DataFrame with added name columns.
        """
        df['Substrate Name'] = df['Substrate InChI'].apply(RpOutput.get_name)
        df['Product Name'] = df['Product InChI'].apply(RpOutput.get_name)
        return df
    
    @staticmethod
    def nodes_to_tsv(df: pd.DataFrame, outdir: str):
        """
        Saves the unique substrate and product names to a TSV file
        Args:
            df (pd.DataFrame): DataFrame containing 'Substrate InChI' and 'Product InChI' columns.
            outdir (str): Directory to save the output TSV file - should end with a '/'
        Returns:
            None
        """
        unique_substrates = pd.Series(df['Substrate InChI'].unique())
        unique_products = pd.Series(df['Product InChI'].unique())

        all_compounds = pd.concat([unique_substrates, unique_products]).unique()
        all_compounds = pd.DataFrame({'inchi': all_compounds})
        all_compounds['name'] = all_compounds['inchi'].apply(RpOutput.get_name)
        all_compounds = all_compounds[['name', 'inchi']]

        all_compounds.to_csv(outdir+'rp_nodes.tsv', sep='\t', index=False, header=True)

    @staticmethod
    def ec_to_name(ec: str) -> str:
        """
        Converts EC number (format "x.x.x.x", not "EC:x.x.x.x") to a human-readable name.
        Queries KEGG to get enzyme names.
        Args:
            ec (str): EC number in the format "x.x.x.x"
        Returns:
            str: Name of the enzyme or 'Unknown' if not found.
        """
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
    
    def split_ec_column(df: pd.DataFrame, ec_column_name: str):
        """
        Split rows with multiple EC numbers (comma-separated) into separate rows.
        
        Args::
            df: pandas DataFrame
            ec_column_name: name of the column containing EC numbers
        
        Returns:
            DataFrame with split rows
        """
        new_rows = []
        
        for idx, row in df.iterrows():
            ec_value = row[ec_column_name]
            
            # Check if EC value is null or doesn't contain comma
            if pd.isna(ec_value) or ',' not in str(ec_value):
                new_rows.append(row)
            else:
                # Convert to string and handle brackets if present
                ec_string = str(ec_value).strip()
                
                # Remove brackets if present (e.g., "[EC1, EC2]" -> "EC1, EC2")
                if ec_string.startswith('[') and ec_string.endswith(']'):
                    ec_string = ec_string[1:-1]
                
                # Split by comma and strip whitespace
                ec_numbers = [ec.strip() for ec in ec_string.split(',')]
                
                # Create a new row for each EC number
                for ec in ec_numbers:
                    new_row = row.copy()
                    new_row[ec_column_name] = ec
                    new_rows.append(new_row)
        
        return pd.DataFrame(new_rows).reset_index(drop=True)
        
    @staticmethod
    def edges_to_tsv(df: pd.DataFrame, outdir: str):
        """
        Saves the edges (reactions) to a TSV file with kegg-like configuration.
        Args:
            df (pd.DataFrame): DataFrame containing the retropath results.
            outdir (str): Directory to save the output TSV file - should end with a '/'
        Returns:
            None
        """
        df['EC number'] = df['EC number'].apply(lambda x: x.strip('[]'))
        df['enzyme_name'] = df['EC number'].apply(lambda x: x if x == 'NOEC' else 'Unknown')
        #df['source_name'] = 'Unknown'
        df['source_name'] = df['Substrate InChI']
        #df['target_name'] = 'Unknown'
        df['target_name'] = df['Product InChI']
        df['color'] = 'RetroPath'
        #df['RPscore'] = df['Score'] + 2  # adding RP penalty for later weighting

        df = df.rename(columns={
            'EC number': 'EC', 
            'Substrate InChI': 'source_inchi', 
            'Product InChI': 'target_inchi', 
            'Substrate SMILES': 'source_smiles',
            'Product SMILES': 'target_smiles',
            'Transformation ID': 'reaction'})
        
        df = RpOutput.split_ec_column(df, "EC")
        
        df.to_csv(outdir+'rp_edges.tsv', sep='\t', index=False, header=True)


    @staticmethod
    def tsv_generator(results_path: str, outdir:str):
        """ 
        Generates edges tsv compatible with kegg edges 
        Args:
            results_path (str): Path to the retropath results.csv file
            outdir (str): Directory to save the output TSV files - should end with a '/'
        Returns:
            None
        """

        results_dir = Path(results_path).parent  # Get the directory containing results.csv
        scope_files = list(results_dir.glob("*scope.csv"))  # Look for any file ending with "scope.csv" in the results directory
        
        if scope_files:
            # Use the first scope file found
            # TODO: add logic to handle multiple matches using source name
            results = scope_files[0]
        else:
            # Fall back to the original results.csv file
            results = results_path

        df = RpOutput.load_res(str(results))
        RpOutput.edges_to_tsv(df, outdir)