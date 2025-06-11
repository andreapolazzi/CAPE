import subprocess
from datetime import datetime
import argparse
import pandas as pd

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


def test_one_to_one(source_name, source_file, sink_files, min_step=1, max_step=5, server=False):
    
    log_filename = f"data_files/log_oto_{source_name}_server.txt" if server else f"data_files/log_oto_{source_name}_laptop.txt"
    with open(log_filename, "w") as log_file:
        log_file.write(f"Log file to test performance on {'server' if server else 'laptop'} \n")
        log_file.write(f"Test executed on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    
    for i in range(min_step, max_step + 1):

        sink_file = sink_files[i-1]
        rules_file = "../retrorules_rr02_rp2_flat_forward.csv"
        outdir = "results/"
        #topx = "100"
        #topx="200"
        topx = "500"
        dmin = "6"
        dmax = "16"

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
            "--msc_timeout", "180"
        ]

        log_message = f"Testing one to one performance for source {source_name} with max_steps={i}\n"
        
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
                df = pd.read_csv("results/results.csv")
                line_count = len(df)
                in_sink_count = df['In Sink'].value_counts().get(1, 0)
                sink_names = df[df['In Sink'] == 1]['Sink name'].unique()
                
                log_file.write(f"Number of lines in results.csv: {line_count}\n")
                log_file.write(f"Number of compounds found in sink: {in_sink_count}\n")
                log_file.write(f"Sink compounds found: {', '.join(sink_names)}\n")
            except FileNotFoundError:
                log_file.write("results.csv file not found\n")


def test_diameters(source_name, source_file, sink_files, min_d=6, max_d=16, min_step=1, max_step=5, server=False):
    
    if args.out_file:
        log_filename = f"data_files/log_diam_{source_name}_server{str(args.out_file)}.txt" if server else f"data_files/log_diam_{source_name}_laptop{str(args.out_file)}.txt"
    else:
        log_filename = f"data_files/log_diam_{source_name}_server.txt" if server else f"data_files/log_diam_{source_name}_laptop.txt"
    with open(log_filename, "w") as log_file:
        log_file.write(f"Log file to test performance on {'server' if server else 'laptop'} \n")
        log_file.write(f"Test executed on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    
    for j in range(min_d, max_d + 1):

        with open(log_filename, "a") as log_file:
            log_file.write("\n" + "="*100 + "\n")
            log_file.write(f"\nTesting diameter from {min_d} to {j}\n")

        for i in range(min_step, max_step + 1):
            with open(log_filename, "a") as log_file:
                log_file.write("\n" + "-"*100 + "\n")
                log_file.write(f"\nTesting max_steps={i}\n")

            sink_file = sink_files[i-1]
            rules_file = "../retrorules_rr02_rp2_flat_forward.csv"
            outdir = "results/"
            topx = "100"
            #topx="200"
            #topx = "500"
            dmin = str(min_d)
            dmax = str(j)

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
                "--msc_timeout", "180"
            ]

            with open(log_filename, "a") as log_file:
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
                    df = pd.read_csv("results/results.csv")
                    line_count = len(df)
                    in_sink_count = df['In Sink'].value_counts().get(1, 0)
                    sink_names = df[df['In Sink'] == 1]['Sink name'].unique()
                    
                    log_file.write(f"Number of lines in results.csv: {line_count}\n")
                    log_file.write(f"Number of compounds found in sink: {in_sink_count}\n")
                    log_file.write(f"Sink compounds found: {', '.join(sink_names)}\n")
                except FileNotFoundError:
                    log_file.write("results.csv file not found\n")
       		 

        
        




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run performance tests for RetroPath2')
    parser.add_argument('-s', '--server', action='store_true', help='Run tests in server mode')
    parser.add_argument('--min_steps', type=int, default=1, help='Minimum number of steps')
    parser.add_argument('--max_steps', type=int, default=5, help='Maximum number of steps')
    parser.add_argument('--ota', action='store_true', help='Run one-to-all test')
    parser.add_argument('--otc', action='store_true', help='Run one-to-central test')
    parser.add_argument('--oto', action='store_true', help='Run one-to-one test')
    parser.add_argument('--dia', action='store_true', help='Run incremental diameters on one-to-one test')
    parser.add_argument('--min_diam', type=int, default=6, help='Minimum diameter for one-to-one test')
    parser.add_argument('--max_diam', type=int, default=16, help='Maximum diameter for one-to-one test')
    parser.add_argument('--out_file', type=int, help='Specify active outputs number')
    #parser.add_argument('--target_ids', type=str, nargs='+', help='Target IDs for one-to-one test')
    args = parser.parse_args()
    
    if args.ota:
        test_one_to_all(min_step=args.min_steps, max_step=args.max_steps, server=args.server)
    
    if args.otc:
        test_one_to_central(min_step=args.min_steps, max_step=args.max_steps, server=args.server)
    
    if args.oto:
        """ if not args.source_id or not args.target_ids:
            print("For one-to-one test, please specify --source_id and --target_ids")
            parser.print_help()
            exit(1) """
        source = "data_files/met.csv"
        name = "MET"
        sinks = [
            "data_files/terephtalate.csv",
            "data_files/ter_dihydrodiol.csv",
            "data_files/dihydroxy_benzoate.csv",
            "data_files/beta_carboxy_muconate.csv",
            "data_files/muconolactone.csv",
            "data_files/dihydrofuran.csv",
            "data_files/oxoadipate.csv",
            "data_files/oxoadipyl_coa.csv",
            "data_files/succinyl_coa.csv"
        ]
        test_one_to_one(source_name=name, source_file=source, sink_files=sinks, min_step=args.min_steps, max_step=args.max_steps, server=args.server)
    
    if args.dia:

        source = "data_files/met.csv"
        name = "MET"
        sinks = [
            "data_files/terephtalate.csv",
            "data_files/ter_dihydrodiol.csv",
            "data_files/dihydroxy_benzoate.csv",
            "data_files/beta_carboxy_muconate.csv",
            "data_files/muconolactone.csv",
            "data_files/dihydrofuran.csv",
            "data_files/oxoadipate.csv",
            "data_files/oxoadipyl_coa.csv",
            "data_files/succinyl_coa.csv"
        ]

        test_diameters(source_name=name, source_file=source, sink_files=sinks, min_d=args.min_diam, max_d=args.max_diam, min_step=args.min_steps, max_step=args.max_steps, server=args.server)


    if not (args.ota or args.otc or args.oto or args.dia):
        print("Please specify either --ota, --otc, --oto or --dia")
        parser.print_help()
        exit(1)
        
    print("Performance tests completed.")
    print(f"Check data_files/log_ot{'a' if args.ota else 'c'}_{('server' if args.server else 'laptop')}.txt for details.")
    print("Check results/ for output files.")
