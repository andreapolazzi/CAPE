import subprocess

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
