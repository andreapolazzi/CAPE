import subprocess
import os
import platform

def launch_example():
    # Determine the operating system
    current_os = platform.system().lower()

    # Define the base directories
    codes_dir = os.path.dirname(os.path.abspath(__file__))  # Directory of the current script
    iGEM_root = os.path.dirname(codes_dir)  # Parent directory of 'codes'
    project_root = os.path.dirname(iGEM_root)  # Parent directory of 'iGEM' repo
    output_dir = os.path.join(iGEM_root, 'results')
    os.makedirs(output_dir, exist_ok=True)  # Create output directory if it doesn't exist

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
    rules_file = os.path.join(retropath_folder, 'example_data', 'rules.csv')
    source_file = os.path.join(retropath_folder, 'example_data', 'source.csv')
    sink_file = os.path.join(retropath_folder, 'example_data', 'sink.csv')
    pathway_length = '4'

    try:
        # Construct the command
        command = command_prefix + [
            retropath_cli, 
            rules_file, 
            source_file, 
            sink_file, 
            pathway_length,
            '--output-dir',
            output_dir
        ]

        # Print paths and OS info for debugging
        print(f"Operating System: {current_os}")
        print("RetroPath Folder:", retropath_folder)
        print("RetroPath CLI Path:", retropath_cli)
        print("Command to be executed:", ' '.join(command))

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