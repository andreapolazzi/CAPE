import subprocess
from datetime import datetime
import argparse

def run_command(command):
    """
    Execute a command in the Linux terminal

    Args: command (list): The command to execute as a list of strings

    Returns: tuple: (return_code, stdout, stderr)
    """
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
    

def test_performance(rp_command, clean=True):
    """
    Test the performance of a RetroPath2 command using the time utility.
    Args:
        rp_command (str): The RetroPath2 command to execute.
    Returns:
        tuple: (return_code, stdout, stderr)
    """
    if clean:
        print("Cleaning up previous results...")
        subprocess.run("rm -r results/*", shell=True)

    command = [
        "/usr/bin/time",
        "-v"
    ]
    command.extend(rp_command)  # Add each element of rp_command to the command list
    print("Running performance test with command:", ' '.join(command))

    return_code, stdout, stderr = run_command(command)
    if return_code != 0:
        print("Error executing command:", stderr)
    else:
        print("Command executed successfully:", stdout)

    return return_code, stdout, stderr



def test_one_to_all(min_step=1, max_step=5, server=False):

    log_filename = "data_files/log_ota_server.txt" if server else "data_files/log_ota_laptop.txt"
    with open(log_filename, "w") as log_file:
        log_file.write(f"Log file to test performance on {'server' if server else 'laptop'} \n")
        log_file.write(f"Test executed on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    sink_file = "data_files/biocyc_cmpd_rpcleaned.csv"
    rules_file = "../retrorules_rr02_rp2_flat_forward.csv"
    outdir = "results/"
    source_file = "data_files/met.csv"
    topx = "100"
    dmin = "6"
    dmax = "16"

    for i in range(min_step, max_step + 1):
        rp_command = [
            "python", "-m", "retropath2_wrapper",
            sink_file,
            rules_file,
            outdir,
            "--source_file", source_file,
            "--max_steps", str(i),
            "--topx", topx,
            "--dmin", dmin,
            "--dmax", dmax,
        ]

        log_message = f"Testing one to all performance for max_steps={i}\n"
        
        with open(log_filename, "a") as log_file:
            log_file.write(log_message)
            return_code, stdout, stderr = test_performance(rp_command)
            log_file.write("STDOUT:\n")
            log_file.write(stdout + "\n")
            log_file.write("STDERR:\n")
            log_file.write(stderr + "\n")
            if return_code:
                log_file.write(f"Error executing command (exit code {return_code})\n")
            else:
                log_file.write("Command executed successfully.\n")
                
            try:
                with open("results/results.csv", "r") as f:
                    line_count = sum(1 for line in f)
                log_file.write(f"Number of lines in results.csv: {line_count}\n")
            except FileNotFoundError:
                log_file.write("results.csv file not found\n")

            log_file.write("\n")
            log_file.write("-" * 80 + "\n")



def test_one_to_central(min_step=1, max_step=5, server=False):

    log_filename = "data_files/log_otc_server.txt" if server else "data_files/log_otc_laptop.txt"
    with open(log_filename, "w") as log_file:
        log_file.write(f"Log file to test performance on {'server' if server else 'laptop'} \n")
        log_file.write(f"Test executed on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    sink_file = "data_files/gem_central_rp.csv"
    rules_file = "../retrorules_rr02_rp2_flat_forward.csv"
    outdir = "results/"
    source_file = "data_files/met.csv"
    topx = "100"
    dmin = "6"
    dmax = "16"
    tout = "30"

    for i in range(min_step, max_step + 1):
        rp_command = [
            "python", "-m", "retropath2_wrapper",
            sink_file,
            rules_file,
            outdir,
            "--source_file", source_file,
            "--max_steps", str(i),
            "--topx", topx,
            "--dmin", dmin,
            "--dmax", dmax,
            "--msc_timeout", tout,
        ]

        log_message = f"Testing one to central metabolites performance for max_steps={i}\n"
        
        with open(log_filename, "a") as log_file:
            log_file.write(log_message)
            return_code, stdout, stderr = test_performance(rp_command)
            log_file.write("STDOUT:\n")
            log_file.write(stdout + "\n")
            log_file.write("STDERR:\n")
            log_file.write(stderr + "\n")
            if return_code:
                log_file.write(f"Error executing command (exit code {return_code})\n")
            else:
                log_file.write("Command executed successfully.\n")
                
            try:
                with open("results/results.csv", "r") as f:
                    line_count = sum(1 for line in f)
                log_file.write(f"Number of lines in results.csv: {line_count}\n")
            except FileNotFoundError:
                log_file.write("results.csv file not found\n")

            log_file.write("\n")
            log_file.write("-" * 80 + "\n")




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run performance tests for RetroPath2')
    parser.add_argument('-s', '--server', action='store_true', help='Run tests in server mode')
    parser.add_argument('--min_steps', type=int, default=1, help='Minimum number of steps')
    parser.add_argument('--max_steps', type=int, default=5, help='Maximum number of steps')
    parser.add_argument('--ota', action='store_true', help='Run one-to-all test')
    parser.add_argument('--otc', action='store_true', help='Run one-to-central test')
    args = parser.parse_args()
    
    if args.ota:
        test_one_to_all(min_step=args.min_steps, max_step=args.max_steps, server=args.server)
    if args.otc:
        test_one_to_central(min_step=args.min_steps, max_step=args.max_steps, server=args.server)
    if not (args.ota or args.otc):
        print("Please specify either --ota or --otc")
        parser.print_help()
        exit(1)
        
    print("Performance tests completed.")
    print(f"Check data_files/log_ot{'a' if args.ota else 'c'}_{('server' if args.server else 'laptop')}.txt for details.")
    print("Check results/ for output files.")