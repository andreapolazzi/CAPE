import subprocess
import os
import platform

def launch_retropath(**kwargs):
    
    # Define the base directories
    codes_dir = os.path.dirname(os.path.abspath(__file__))  # Directory of the current script
    iGEM_root = os.path.dirname(codes_dir)  # Parent directory of 'codes'
    project_root = os.path.dirname(iGEM_root)  # Parent directory of 'iGEM' repo
    if "output_dir" in kwargs:
        output_dir = kwargs["output_dir"]  # output directory if provided by user
    else:
        output_dir = os.path.join(iGEM_root, 'results')
    os.makedirs(output_dir, exist_ok=True)  # Create output directory if it doesn't exist

    # Determine the operating system
    current_os = platform.system().lower()
    # Determine the correct folder and executable based on the operating system
    if current_os == 'windows':
        retropath_folder = os.path.join(project_root, 'retropath_1.0.0-alpha_windows-x64')
        retropath_cli = os.path.join(retropath_folder, 'RetroPath.Cli.exe')
        command_prefix = []  # No prefix needed for Windows executable
    elif current_os == 'darwin':  # macOS
        retropath_folder = os.path.join(project_root, 'osx-arm64')
        retropath_cli = os.path.join(retropath_folder, 'RetroPath.Cli.dll')
        command_prefix = ['dotnet']  # Requires dotnet to run .dll
    elif current_os == 'linux':
        retropath_folder = os.path.join(project_root, 'retropath_1.0.0-alpha_linux-x64')
        retropath_cli = os.path.join(retropath_folder, 'RetroPath.Cli')
        command_prefix = []  # No prefix needed for Linux
    else:
        print(f"Unsupported operating system: {current_os}")
        return

    # Define full paths to files
    if "rules_file" in kwargs:
        rules_file = kwargs["rules_file"]
    else:
        rules_file = os.path.join(iGEM_root, 'data_files', 'rules.csv')
    if "source_file" in kwargs:
        source_file = kwargs["source_file"]
    else:
        source_file = os.path.join(iGEM_root, 'data_files', 'source.csv')
    if "sink_file" in kwargs:
        sink_file = kwargs["sink_file"]
    else:
        sink_file = os.path.join(iGEM_root, 'data_files', 'sink.csv')
    if "pathway_length" in kwargs:
        pathway_length = kwargs["pathway_length"]
    else:
        pathway_length = '4'
    if "min_diameter" in kwargs:
        min_diameter = kwargs["min_diameter"]
    else:
        min_diameter = '0'
    if "max_diameter" in kwargs:
        max_diameter = kwargs["max_diameter"]
    else:
        max_diameter = '1000'
    if "max_structures" in kwargs:
        max_structures = kwargs["max_structures"]
    else:
        max_structures = '100'
        # --max-structures 

    try:
        # Construct the command
        command = command_prefix + [
            retropath_cli, 
            rules_file, 
            source_file, 
            sink_file, 
            pathway_length,
            '--output-dir',
            output_dir,
            '--min-diameter',
            min_diameter,
            '--max-diameter',
            max_diameter,
            '--max-structures',
            max_structures
        ]

        # Print paths and OS info for debugging
        print(f"Operating System: {current_os}")
        print("RetroPath Folder:", retropath_folder)
        print("RetroPath CLI Path:", retropath_cli)
        print("Command to be executed:", ' '.join(command))

        # Run the command
        print("Launching RetroPath...")
        result = subprocess.run(command, check=True, capture_output=True, text=True)

        # Print the output
        print("Standard Output:")
        print(result.stdout)
        print("Standard Error:")
        print(result.stderr)

    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running RetroPath:")
        print(f"Return Code: {e.returncode}")
        print("Standard Output Process Error:")
        print(e.stdout)
        print("Standard Error Process Error:")
        print(e.stderr)
    except FileNotFoundError:
        print("Error: RetroPath CLI executable not found. Please check the file path.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    launch_retropath(pathway_length='3')