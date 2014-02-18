alignment_filtering
===================

This repository contains materials used for the manuscript, "Limited utility of residue masking for positive-selection inference" by Stephanie J. Spielman, Eric T. Dawson, Claus O. Wilke.

Description of directories and contents as follows.

Manuscript/
	Contains Latex manuscript and corresponding materials (figures, citations, style files, etc)

OSX_scoringprograms/
	Contains source code for C++ programs to carry out residue scoring from bootstrap replicates. All executables have been compiled for Mac OSX. Additionally contains an initial implementation written in python (scoringtest.py) which is no longer used.

SJSGuidance/
	Contains all code and a README file to run our re-implemented Guidance. Note that we call it SJSGuidance due to the author's initials, SJS, and to distintuish from the original Guidance written by Penn et al. (2010). See the README within that directory for details.

data/
	Contains final data files (data/parsed_data/), figures (data/Rcode_plots/), statistics (data/stats) for manuscript

scripts/
	Contains code for analyses conducted in manuscript, including sequence simulation, selection inference with FUBAR and PAML M8, and parsing of the raw data files. Note that, as alignments were conducted with our Guidance reimplementation, you should see the directory SJSGuidance for alignment code.