import subprocess

class RP():

    @staticmethod
    def rp_source_file(sink: str, rules:str, outdir:str,
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
            return 1, "", str(e)
        
    @staticmethod
    def rp_single_source(sink: str, rules:str, outdir:str, 
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
            return 1, "", str(e)