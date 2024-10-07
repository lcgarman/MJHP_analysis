# MJHP_analysis
The Mott-Jones Hamilton Population analysis. A program for analyzing the Mott-Jones effect in intermetallic compounds using abinit.

## Installation
Provided you have an SSH key to clone the repo, run these commands:
```
git@github.com:lcgarman/MJHP_analysis.git
cd MJHP_analysis/src
make -f makefile-mjhpHKL
make -f makefile-mjhp2theta
make -f makefile-prep2theta
```
This will clone the repo into your current directory, change to the necessary directory, and run the makefiles to compile the necessary executables. You may copy the binaries into any directly in your $PATH (/usr/local/bin, for instance). 

## Dependencies

### Compilation

The GNU Scientific Libary, libxc, and fftw3 are required for compiling the MJHP_analysis.

### Running Calculation

The MJHP_analysis requires outputs from abinit version 7.10.5. This method requires HGH (Hartwigsen-Goedecker-Hutter) norm-conserving pseudopotentials and LDA exchange-correlation functionals. The MJHP_analysis expects the necessary abinit output files to be named: FILENAME_o_XXX, where XXX refers to the different file types. I recommend working in two separate directories for the MJHP_analysis, let’s call these directories:

  2theta
  
  HKL
  
Into the 2theta directory copy over the abinit FILENAME.in, FILENAME.files (containing the paths to the pseudopotentials), and FILENAME_i_DEN (the output density from the single-point energy FILENAME_o_DEN will be used as the starting density FILENAME_i_DEN for the non-self consistent calculation to follow). 
Into the HKL directory copy the FILENAME.out, FILENAME_o_POT, and FILENAME_o_WFK files from the single-point energy calculation. 
One additional input file is required in both of these directories: FILENAME.mjin. Start by making this file in the 2theta directory with the following format:

IN_FILENAME

OUT_COMPOUND

vec #

xY

The first line, IN_FILENAME, should give the file string corresponding to your abinit output files. The OUT_COMPOUND is what you would like the starting file string to be of the output files from the MJHP_analysis programs. The third line expects a string followed by a number (double or integer) that corresponds to the number of valence electrons in the unit cell of the compound. The fourth line, xY, should be two letters where the first (x) corresponding to the crystal family and the second to the centering type (Y). This line is case sensitive, so x should be lowercase and Y uppercase, identical to the letters of a Pearson symbol. 
Once the FILENAME.mjin file and the necessary abinit files are within the 2theta directory run the executable:
```
prepare_mjhp2theta
```
This will find the high symmetry k-points of interest and copy them into your abinit FILENAME.in file, as well as set up an OUT_COMPOUND.rflc file necessary for the MJHP_analysis to follow. Run the non-self consistent ABINIT calculation with this modified FILENAME.in. Following the completion of this job, in the same directory run this executable to perform the first MJHP_analysis:
```
mjhp2theta_analyze.
```
This will generate one OUT_FILENAME.mjout file that can be plotted using matlab program:
```
plot_MJHP_2theta(OUT_COMPOUND_2theta.mjout)
```
This completes the first MJHP analysis. 

Now, navigate to the HKL directory, and copy over the FILENAME.mjin file from the 2theta directory. This file now requires a few additional lines as follows:

IN_FILENAME

OUT_COMPOUND

vec #

xY

no_HKL #

HKL

h1 k1 l1

h2 k2 l2

hn kn ln

The first four lines should not change from the 2theta to HKL calculations. The fifth line, should be a string, followed by an integer of how many datasets you want to run. The sixth line should be a string “HKL”. Following this line should be a list of the conventional HKL indices you are interested in analyzing. The number of lines following HKL should correspond to the integer following the no_HKL on the fifth line. Following this set up, run the executable:
```
mjhpHKL_analyze
```
This will generate OUT_COMPOUND_h1k1l1.mjout (HKL corresponding to the conventional indices of the reflections of interest) files for each dataset will be generate. The results from this analysis may be plotted using 
```
plot_MJHP_HKL(OUT_COMPOUND_h1k1l1.mjout, minY, maxY, directory)
```
The  minY and maxY correspond to the min and max energy range you want to view, and directory can either be a path to where to save the figure files (.tif and .fig), or 0 to not save the figure files. 

This completes the MJHP_analysis. Below is a sample mjin file for Cu5Zn8:

Cu5Zn8

Cu5Zn8_gammabrass

vec 42.0

cI

no_HKL 3

HKL

3 3 0


4 1 1 

-4 1 1

## Citation 
If you use the MJHP_analysis in your research, please cite the following paper:

(1) Garman, L. C.; Fredrickson, D. C. The Mott-Jones Hamilton Population: Energetics of Fermi Sphere–Jones Zone Interactions in Intermetallic Structures. The Journal of Physical Chemistry C 2024, 128 (34), 14442-14457. DOI: 10.1021/acs.jpcc.4c03450.

