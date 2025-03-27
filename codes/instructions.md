## expected folder organization:

```
project_root/
│
├── codes/
│   ├── prova1.py
│   └── prova1_crossplat.py 
|
└── retropath_1.0.0-alpha_linux-x64/ (for Linux)
    ├── RetroPath.Cli (or .exe or .dll depending on OS)
    └── example_data/
        ├── rules.csv
        ├── source.csv
        └── sink.csv
```

The folder for retropath changes based on OS:  
    - retropath_1.0.0-alpha_windows-x64.zip/ (for Windows) <br>
    - retropath_1.0.0-alpha_osx-arm64.zip/ (for macOS)

codes/ directory will contain all codes needed for the project. <br>
project_root/ is the github repo iGEM. <br>
To use retropath but not only.

## .gitignore

Create a .gitignore in the folder codes/ 
(see https://www.w3schools.com/git/git_ignore.asp for more info): <br>
1. execute the command: touch .gitignore  
2. open it and write inside:
```
# ignore results folder, as it changes each time  scripts are executed

results/
```
This is needed to avoid uploading the results of each run to github 