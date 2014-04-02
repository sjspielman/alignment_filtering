/* Automates FUBAR with the following input. Written 3/12/13 SJS.
MUST have FUBAR.bf and the directory FUBAR_HBL in working directory. */

BASEDIR = "/state/partition1/sjs3495/fubarCPU-266314/FUBARmaterials/";
//BASEDIR="/home/sjs3495/FUBARmaterials/";
datafile="temp.fasta";
treefile="tree.tre";

inputRedirect = {};
inputRedirect["01"]="Universal";         //Genetic code
inputRedirect["02"]="1";                 //Number of files to process
inputRedirect["03"]=BASEDIR+datafile;    //Fasta file, full path
inputRedirect["04"]=BASEDIR+treefile;    //Tree file, full path
inputRedirect["05"]="20";                //GrfubarCPU-266314 size, DxD grfubarCPU-266314. Default=20
inputRedirect["06"]="5";                 //Number MCMC chains, default=5
inputRedirect["07"]="2000000";            //Length of each MCMC chain. default is 2mil
inputRedirect["08"]="1000000";            //Burnin for each MCMC chain. default 1 mil
inputRedirect["09"]="100";               //Number of samples to draw
inputRedirect["10"]="0.5";               //Dirichlet prior parameter, this is the default.

ExecuteAFile ("FUBAR.bf", inputRedirect);
