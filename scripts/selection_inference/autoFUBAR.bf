/* Automates FUBAR with the following input. Written 4/18/14 SJS.
Note that this version considers only dN variation (dS set to 1).
MUST have FUBAR.bf and the directory FUBAR_HBL in working directory (BASEDIR). */

BASEDIR = "placeholder/";
datafile="temp.fasta";
treefile="tree.tre";

inputRedirect = {};
inputRedirect["01"]="Universal";         //Genetic code
inputRedirect["02"]="1";                 //Number of files to process
inputRedirect["03"]=BASEDIR+datafile;    //Fasta file, full path
inputRedirect["04"]=BASEDIR+treefile;    //Tree file, full path
inputRedirect["05"]="No";                  //Consider only variation in dN
inputRedirect["06"]="100";                //100 grid points along the 1D dN axis
inputRedirect["07"]="5";                 //Number MCMC chains, default=5
inputRedirect["08"]="2000000";            //Length of each MCMC chain. default is 2mil
inputRedirect["09"]="1000000";            //Burnin for each MCMC chain. default 1 mil
inputRedirect["10"]="100";               //Number of samples to draw
inputRedirect["11"]="0.5";               //Dirichlet prior parameter, this is the default.

ExecuteAFile ("FUBAR.bf", inputRedirect);
