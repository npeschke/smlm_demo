# Demo repository for the SMLM package

## Installation Guide for Linux
1. open a terminal 
2. create a conda environment
   - `conda create -n voronoi_analysis python=3.9`
3. activate the environment
   - `conda activate voronoi_analysis`
4. install the package via pip
   - `pip install git+https://github.com/npeschke/smlm`

The typical install time is less than 10 min on a "normal" desktop computer

## Running the demo
1. clone the repository to your machine
   - `git clone https://github.com/npeschke/smlm_demo.git`
2. open the repo
   - `cd smlm_demo`
3. run the demo
   - `python demo.py`

This will read the localization files from the `data/orte` directory,
run the analysis and plot figures similar to the ones in Figure 5
in the paper that are then located in the `results` directory.
The demo should run in under 10 minutes on a "normal" desktop
computer.

## Instructions for use
The demo repository also contains the notebooks used to create
the plots in the paper. Due to filesize restrictions in github the
data needed to reproduce them is not included but available upon
reasonable request. For their use, `jupyter` needs to be installed
in the environment with the following command.
```bash
conda install jupyter
```
