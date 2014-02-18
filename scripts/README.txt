SJS 2/17/14. Scripts used to perform analysis in "Limited utility of residue masking for positive selection inference." authors SJS, ETD, COW.

parsing/
	Used to parse output files from fubar and paml selection inferences
		parsefubar.py and parsepaml.py used to collect positively selected sites at a posterior probability cutoff of 0.9. Generates TPR, FPR, etc
		sweepfubar.py and sweeppaml.py used to cycle over posterior probability cutoffs 0-100 and generate TPR, FPR, etc. Useful for ROC analysis and relationship between TPR/FPR and posterior probability.

selection_inference/
	Basic scripts for conducting PAML (codeml.ctl) and FUBAR (autoFUBAR.bf) inferences

simulation/
	Scripts used to simulate sequences with Indelible: runindelible.py and control.txt
	trees/ 
		See manuscript for phylogeny citations.
			or5: 11 sequences
			rho: 26 sequences
			prk: 60 sequences
			flat:158 sequences
	HA/ 
		Inference of H1N1 HA parameters