# MJHP_analysis
The Mott-Jones Hamilton Population analysis. A program for analyzing the Mott-Jones effect in intermetallic compounds using ABINIT.

## Installation
Provided you have an SSH key to clone the repo, run these commands:
```
cd MJHP_analysis/src
make -f makefile-mjhpHKL
make -f makefile-mjhp2theta
make -f makefile-prep2theta
```
This will clone the repo into your current directory, change to it, and run the makefiles to compile the necessary executables. You may copy the binaries into any directly in your $PATH (/usr/local/bin, for instance). 

## Dependencies
### Compilation
The GNU Scientific Libary, libxc, and fftw3 are required for compiling the MJHP_analysis.
### Running Calculation
The MJHP_analysis requires outputs from abinit version 7.10.5. This method requires HGH (Hartwigsen-Goedecker-Hutter) norm-conserving pseudopotentials, and LDA exchange-correltation functionals. The MJHP_analysis expects the ABINIT output files to be named: FILENAME_o_XXX, where XXX refers to different file types. I recommend working in two separate directories following the completion of the ABINIT single-point energy calculation, lets call these directories:
  2theta
  HKL
Into the 2theta directory copy over the abinit FILENAME.in, FILENAME.out, and FILENAME_i_DEN (the output density from the single-point energy FILENAME_o_DEN will be used as the starting density FILENAME_i_DEN for the non-self consitent calculation to follow). Into the HKL directory copy over the FILENAME.out, FILENAME_o_POT, FILENAME_o_WFK files. One additional input file is required in both directories: FILENAME.mjin. Start by making this file in the 2theta directory with the following format:
IN_FILENAME
OUT_COMPOUND
vec #
xY

The first line, IN_FILENAME, should give the string for your abinit output files. The OUT_COMPOUND is what you would like the starting filestring to be of the output files from the MJHP_analysis programs. The third line expects a string followed by a number (double or integer) that corresponds to the number of valence electrons in the unit cell of the compound. The fourth line, xY, should be two letters the first (x) corresponding to the crystal family and the second to the centering type (Y). This line is case sensitive, so x should be lowercase and Y uppercase, identical to the letters of a Pearson symbol. 

Once the FILENAME.mjin file and the necessary abinit files are within the 2theta directory run the executable:
prepare_mjhp2theta
This will find the high symmetry k-points of interest and copy them into your abinit FILENAME.in file, as well as set up an OUT_COMPOUND.rflc file necessary for the MJHP_analysis to follow. Run the non-self consistent ABINIT calculation with this modified FILENAME.in file. Following the completion of this job, in the same directory run the following executable to perform the first MJHP_analysis;
mjhp2theta_analyze.

This will generate one OUT_FILENAME.mjout file that can be plotted using matlab program:
plot_MJHP_2theta.m

This completes the first MJHP analysis. 

Now, navigate the HKL directory, and copy over the FILENAME.mjin file from the 2theta directory to the HKL directory. This file now requires a few additional lines as follows:

IN_FILENAME
OUT_COMPOUND
vec #
xY
no_HKL #
HKL
h1 k1 l1
h2 k2 l2
hn kn ln

The first four lines should not change from the 2theta to HKL calculations. The fifth line, should be a string, "no_HKL" followed by an integer of how many datasets you want to run. The sixth line should be a string HKL. Folloiwng this line should be the CONVENTION HKL indices you are interested in. The number of lines following HKL should be the same as the integer following no_HKL. Following this set up, run the executable:
mjhpHKL_analyze
This will generate OUT_COMPOUND_h1k1l1.mjout (HKL corresponding to the conventional indices of the reflections of interest) files for each dataset will be generate. The results from this analysis may be plotted using 

plot_MJHP_HKL(OUT_COMPOUND_h1k1l1.mjout, minY, maxY, directory). - minY and maxY correspond to the min and max energy range you want to view, and directory can either be a path to where you want to save the figure files (.tif and .fig), or 0 to not save the figure files. 

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





