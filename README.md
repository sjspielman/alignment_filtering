alignment_filtering
===================

This repository contains materials used for the manuscript, "Limited utility of residue masking for positive-selection inference" by Stephanie J. Spielman, Eric T. Dawson, Claus O. Wilke.
Please direct all questions to Stephanie, at stephanie.spielman@gmail.com (or stephanie.spielman@utexas.edu).

Description of directories and contents as follows.

Manuscript/
	Contains Latex manuscript and corresponding materials (figures, citations, style files, etc)

OSX_scoringprograms/
	Contains source code for C++ programs to carry out residue scoring from bootstrap replicates. All executables have been compiled for Mac OSX. Additionally contains an initial implementation written in python (scoringtest.py) which was NOT used in the study.

PhyloGuidance/
	Contains all code and a README file to run our re-implemented Guidance, PhyloGuidance. See the README within that directory for details. Note that the name "PhyloGuidance" doesn't carry over into the manuscript, but is a nicer name for this directory than "Guidance_Reimplementation."

data/
	Contains final parsed data files and statistics. See README within for details.	

scripts/
	Contains code for analyses conducted in manuscript, including sequence simulation, selection inference with FUBAR and PAML M8, and parsing of the raw data files. All alignments were conducted within PhyloGuidance. The README in the scripts/ directory contains more details.