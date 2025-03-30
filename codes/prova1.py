import subprocess
import os
# ciao
def launch_example():
    # Define the base directories
    codes_dir = os.path.dirname(os.path.abspath(__file__))  # Directory of the current script
    project_root = os.path.dirname(codes_dir)  # Parent directory of 'codes'
    retropath_folder = os.path.join(project_root, 'retropath_1.0.0-alpha_linux-x64')

    # Define full paths to files
    retropath_cli = os.path.join(retropath_folder, 'RetroPath.Cli')
    rules_file = os.path.join(retropath_folder, 'example_data', 'rules.csv')
    source_file = os.path.join(retropath_folder, 'example_data', 'source.csv')
    sink_file = os.path.join(retropath_folder, 'example_data', 'sink.csv')
    pathway_length = '4'

    try:
        # Construct the command as it would be in the terminal
        command = [
            retropath_cli, 
            rules_file, 
            source_file, 
            sink_file, 
            pathway_length
        ]

        # Print paths for debugging
        print("RetroPath CLI Path:", retropath_cli)
        print("Rules File Path:", rules_file)
        print("Source File Path:", source_file)
        print("Sink File Path:", sink_file)

        # Run the command
        print("Launching RetroPath with example data...")
        result = subprocess.run(command, check=True, capture_output=True, text=True)

        # Print the output
        print("Standard Output:")
        print(result.stdout)
        
        print("Standard Error:")
        print(result.stderr)

    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running RetroPath:")
        print(f"Return Code: {e.returncode}")
        print("Standard Output:")
        print(e.stdout)
        print("Standard Error:")
        print(e.stderr)
    except FileNotFoundError:
        print("Error: RetroPath CLI executable not found. Please check the file path.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    launch_example()