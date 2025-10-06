# Multi-level 1D TDSE solver
There are several techniques to solve the time-dependent Schrödinger equation (TDSE). Here, we have a one-dimensional TDSE solver employing the split-operator method, one of the efficient spectral methods. The split-operator method is made efficient due to the accessibility of the fast Fourier transform algorithms. In this package, we use Fast Fourier Transform in the west (FFTW3) package. Below, we discuss the general installation and capabilities of our package.

**Note:** A development phase GUI is now available. More on it below.

## Installation
### System requirements and dependencies
This package has only been tested on Linux systems (Ubuntu/Debian). The distribution for Mac and Windows is being worked on.

1. Nix package manager
2. git (optional)

### Quick run on Systems with Nix
If you have a working Nix installation with nix flakes turned on (see below), then you may use the following commands to quickly build and run the package.
```
$ git clone https://github.com/saurabbh14/Multilevel_1D-TDSE.git   # Cloning the GitHub repo
$ cd Multilevel_1D-TDSE                                            # Package directory
$ nix build                                                        # Package build
$ ./results/bin/ML-TDSE input.ini                                  # Run with test input parameters (input.ini)
```
The further usage and test output are given in USAGE.md

### Installing Nix
To install Nix, follow the instructions on the [Nix download](https://nixos.org/download/) page.

**Note:** The Multi-user installation (assuming you have root access) is recommended; however, the single-user installation should work too!

After successfully installing Nix, add the following line to the file /etc/nix/nix.conf.
```
experimental-features = nix-command flakes
```

### Bulding package
After installing Nix, just clone this [github repo](https://github.com/saurabbh14/Multilevel_1D-TDSE.git) to the intended directory.
```
$ git clone https://github.com/saurabbh14/Multilevel_1D-TDSE.git 
```
Or you can simply download the [ZIP file of the repo](https://github.com/saurabbh14/Multilevel_1D-TDSE/archive/refs/heads/master.zip) and unzip it in the intended directory.

Navigate to the directory where you have cloned or unzipped the git repo.
```
$ cd $Dir_path
```

Usually `$Dir_path = Multilevel_1D-TDSE`. Once in the directory, just run the following command
```
$ nix build
```

This will create an executable file at $Dir_path/result/bin/ML-TDSE.

### Test run
To run the test calculations, you can use the provided input file for the H2+ molecular ion.
```
$ ./results/bin/ML-TDSE input.ini
```  

For more in depth documentation of the package is given in [`DOCUMENTATION.md`](./documentation/DOCUMENTATION.md) file and for some usage examples are listed in [`USAGE.md`](./documentation/USAGE.md) file. 

## Graphical User Interface (GUI)

A simple graphical interface is provided for editing input files, running simulations, and visualizing outputs.

**Features:**
- Edit and save `input.ini` with a built-in text editor.
- Change common parameters using "Quick fields".
- Build the project (`nix build`) and run simulations with one click.
- View simulation output and logs in real time.
- Plot 1D results (e.g. `density_1d.out`) directly in the GUI.

**Requirements:**  
- Python 3.8+  
- PySide6, matplotlib, numpy

**Install dependencies:**
Creating a vitual python environment in the installation directory (`$Dir_path`) is recommended. You may use the following commands to do so or check this [venv documentation](https://docs.python.org/3/library/venv.html): 
```sh
# Assuming we are in $Dir_path
python3 -m venv ".venv"  # will create .venv directory in $Dir_path directory
```
And then install the dependencies in the virtual environment:
```sh
source .venv/bin/activate  # in Linux
python3 -m pipx install -r tools/requirements.txt
```

**Launching GUI**
Currently, there are two flavours of GUIs. We recommend to use the gui based in `tools/pyqt_gui.py`. And it can be launched from anywhere using following command:
```sh
python3 $Dir_path/tools/pyqt_gui.py
```

For more information about the GUI user guide check [`GUI_DOC.md`](./documentation/GUI_DOC.md)

