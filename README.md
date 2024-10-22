# Multi-level 1D TDSE solver
There are several techniques to solve time-dependent Schrödinger equation (TDSE). Here, we have a one-dimensional TDSE solver employing split-operator method, one of the efficient spectral methods. The split-operator method is made efficient due the accessiblity of the fast fourier transform algorithms. In this package we use Fast Fourier Transform in the west (FFTW3) package. Below we discuss the general installation and capabilities of our package.
 
## Installation
### System requirements and dependencies
This package has only been tested on Linux systems (Ubuntu/debian). The distribution for Mac and Windows is being worked on.

1. Nix package manager
2. git (optional)

### Installing Nix
To install Nix, follow the instructions on the [Nix download](https://nixos.org/download/) page.

**Note:** The Multi-user installation (assuming you have the root access) is recommended, however the single-user installation should work too!

After successfully installing Nix, add following line to the file /etc/nix/nix.conf.
```
experimental-features = nix-command flakes
```

### Bulding package
After installing Nix, just clone this [github repo](https://github.com/saurabbh14/Multilevel_1D-TDSE.git) to the intended directory.
```
$ git clone https://github.com/saurabbh14/Multilevel_1D-TDSE.git 
```
Or you can simply downlod the [ZIP file of the repo](https://github.com/saurabbh14/Multilevel_1D-TDSE/archive/refs/heads/master.zip) and unzip in the intended directory.

Navigate to the directory where you have cloned or unzip the git repo.
```
$ cd $Dir_path
```

Once in the directory, just run following command
```
$ nix build
```

This will create an executable file at $Dir_path/result/bin/ML-TDSE.

### Test run
To run, the test calcuations, you can use the provided input file for H2+ molecular ion.
```
./results/bin/ML-TDSE input.ini
```  
 
