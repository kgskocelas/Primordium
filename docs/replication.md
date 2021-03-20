# Replicating experiments

Here we'll walk you through how to replicate the expriments for the paper. 

These instructions were designed for a Unix-based system. Windows compilation is possible, but would require significant reworking. 

## Download
All you need to replicate this project is included in this GitHub repository. 
First, make sure you have [git installed](https://github.com/git-guides/install-git).

Once git is installed on your machine, navigate to your preferred location and clone the repo. 
Via ssh: 
```
git clone git@github.com:kgskocelas/SpatialRestraint.git <DIR_NAME> --recurse-submodules
```
Via https:
```
git clone https://github.com/kgskocelas/SpatialRestraint.git <DIR_NAME> --recurse-submodules
```
Where you replace <DIR_NAME> with whatever you want the directory to be named, or delete it to name the directory after the repo.

## Prerequisites
This project requires serveral dependencies to compiler and run. 
The following subsections will walk you through the installation process for each of them. 
NOTE: If you used the `--recure-submodules` flag when you cloned the repo, you can skip to the R Modules section!

### Empirical
The C++ core of this project comes from the [Empirical library for scientific software](https://github.com/devosoft/Empirical). 
Empirical is setup as a git submodule inside the applications/source directory. 
To download Empirical, you run: 
```
git submodule init
git submodule update
``` 
Please note will also download other repos that we'll cover shortly. 

To download Empirical separately, clone the repo onto your machine in the location of you choice

Via ssh:
```
git clone git@github.com:devosoft/Empirical.git
```
Via https:
```
git clone https://github.com/devosoft/Empirical.git
```

Since Empirical is a header-only C++ library, it's compilation will be handled when we compile this repo. 

### Python Modules
Most of the connector scripts are written in Python 3. 
The only external dependency is pyvarco, which makes it easier to do combinatorial expansions of parameter sets. 

Pyvarco would have been included in the git submodules if you downloaded Empirical that way. 

If you'd rather install pyvarco by hand, simply run the following: 
```
python3 -m pip install pyvarco
```

### R Libraries
All of our data anlayses are conducted in R, and unfortunately we cannot easily include the required R libraries. 

To install the required libraries, open R in your target environment and run the following lines: 
```{r}
install.packages('ggplot2')
install.packages('dplyr')
install.packages('hash')
install.packages('khroma')
install.packages('scales')
install.packages('ggridges')
```

## Compilation
With all the external software in place, you are ready to compile! 

To start, navigate to the `application` directory
```{bash}
cd application
```
If you downloaded Empirical manually, you will need to edit the Makefile in this directory to point to the include directory of your Empirical install. The important line is right at the top, marked by `EMP_DIR=`

Once that is done (or you downloaded Empirical via git submodules) simple run the Makefile to compile: 
```{bash}
make
```
This is assuming you have make installed, as well as g++ 9.0+.

If you compiled without errors, you can now run the app with
```{bash}
./bin/SpatialRestraint
```


