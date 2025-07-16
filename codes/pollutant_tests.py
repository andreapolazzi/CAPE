import subprocess
import pandas as pd
import time

def rp_launch(sink_file, rules_file, outdir, source_file, steps=3, topx=100, dmin=6, dmax=16):
    
    command = [
            "python", "-m", "retropath2_wrapper",
            sink_file,
            rules_file,
            outdir,
            "--source_file", source_file,
            "--max_steps", str(steps),
            "--topx", topx,
            "--dmin", dmin,
            "--dmax", dmax,
        ]
    
    try:
        # Print command for debugging
        print("Executing command:", ' '.join(command))
        
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
    
def rp_single_source(sink_file, rules_file, outdir, source_name, source_inchi, steps=3, topx=100, dmin=6, dmax=16):
    
    command = [
            "python", "-m", "retropath2_wrapper",
            sink_file,
            rules_file,
            outdir,
            "--source_name", source_name,
            "--source_inchi", source_inchi,
            "--max_steps", str(steps),
            "--topx", str(topx),
            "--dmin", str(dmin),
            "--dmax", str(dmax)
        ]
    
    try:
        # Print command for debugging
        print("Executing command:", ' '.join(command))
        
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
    

def activate_rp_env():
    """ not been tested """
    result = subprocess.run(
        ["conda", "activate", "retropath_env"],
        capture_output=True,
        text=True,
        shell=True  # Use shell=True to activate conda environment
    )
    return result.returncode, result.stdout, result.stderr

def deactivate_rp_env():
    """ not been tested """
    result = subprocess.run(
        ["conda", "deactivate"],
        capture_output=True,
        text=True,
        shell=True  # Use shell=True to deactivate conda environment
    )
    return result.returncode, result.stdout, result.stderr


if __name__ == "__main__":

    df = pd.read_csv("./data_files/pollutants.csv")

    sink_file = "data_files/biocyc_cmpd_rpcleaned.csv"
    rules_file = "../retrorules_rr02_rp2_flat_forward.csv"

    session_num = input("Enter session number: ")
    outdir = f"results/screen_{session_num}"
    
    source_name = df.iloc[int(session_num)-1]['Name']
    source_inchi = df.iloc[int(session_num)-1]['InChI']

    steps = int(input("Enter number of steps (default 3): ") or 3)

    returncode, stdout, stderr = rp_single_source(sink_file, rules_file, outdir, source_name, source_inchi, steps=steps)
    
    print("Return code:", returncode)
    print("Standard Output:", stdout)
    print("Standard Error:", stderr)