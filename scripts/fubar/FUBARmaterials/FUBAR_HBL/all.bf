fprintf (stdout, "[FUBAR PHASE 1] Optimizing relative branch lengths under the nucleotide REV model\n");

vectorOfFrequencies = overallFrequencies;
// 'overallFrequencies' is supplied by _MFReader_.ibf

SKIP_HARVEST_FREQ = 1;
LoadFunctionLibrary ("GRM", {"00":"Global"});

populateTrees ("nuc_tree", fileCount);
ExecuteCommands(constructLF ("nucLF", "nucData", "nuc_tree", fileCount));
Optimize (nuc_res, nucLF);

fprintf (stdout, "[FUBAR PHASE 1 FINISHED] log(L) = ", nuc_res[1][0], "\n");
for (k = 1; k <= fileCount; k += 1)
{
    fprintf (stdout, "\tLength of tree ", k, " (substitutions/site) = ", +Eval("BranchLength(nuc_tree_" + k + ",-1)"),"\n");
}

LF_NEXUS_EXPORT_EXTRA = "positionalFrequencies = " + positionFrequencies + ";";
SetDialogPrompt ("Write nucleotide model fit to");
LIKELIHOOD_FUNCTION_OUTPUT = 7;
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, nucLF);
DeleteObject (nucLF);
fscanf  			(stdin,"String", nucFit);
fscanf  			(stdin,"String", codonFit);
fscanf  			(stdin,"String", gridInfoFile);
fscanf              (stdin,"Number", grid_points);

ExecuteAFile        (PATH_TO_CURRENT_BF + "FUBAR_tools.ibf");
prefix              = "existingFit";
timer               = Time (1);
ExecuteAFile        (nucFit,"","existingFit");
ExecuteAFile        (Join(DIRECTORY_SEPARATOR,{{HYPHY_LIB_DIRECTORY[0][Abs(HYPHY_LIB_DIRECTORY)-2],"TemplateBatchFiles","Utility","LocalMGREV.bf"}}));
ExecuteAFile        (Join(DIRECTORY_SEPARATOR,{{HYPHY_LIB_DIRECTORY[0][Abs(HYPHY_LIB_DIRECTORY)-2],"TemplateBatchFiles","TemplateModels","CF3x4.bf"}}));

nuc3x4 = CF3x4 (existingFit.positionalFrequencies,GeneticCodeExclusions);
PopulateModelMatrix ("MGLocalQ", nuc3x4);
vectorOfFrequencies = BuildCodonFrequencies(nuc3x4);
Model MGLocal = (MGLocalQ, vectorOfFrequencies,0);

GetString (nucLF_Info, existingFit.nucLF, -1);
fileCount = Columns (nucLF_Info["Trees"]);

AC := existingFit.AC;
AT := existingFit.AT;
CG := existingFit.CG;
CT := existingFit.CT;
GT := existingFit.GT;

GetString           (mgBrLen, MGLocal, -1);
ExecuteCommands     ("GetString (nucBrLen, "+(nucLF_Info["Models"])[0]+",-1);");
ExecuteCommands     ("GetString (nucMdlInfo, "+(nucLF_Info["Models"])[0]+",-2);");
ExecuteCommands     ("GetString (nuclParamInfo, "+nucMdlInfo["RATE_MATRIX"]+",-1)");

assert              (Columns(nuclParamInfo["Local"])==1,"The nucleotide model must have exactly one local parameter");
paramName       =   (nuclParamInfo["Local"])[0]; 
ExecuteCommands     ("FindRoot(nf,`nucBrLen`-1,`paramName`,0,1e10);");
paramName       =   paramName[Abs(prefix) + 1][Abs(paramName)-1];

//----------------------------------------------------------------------------


rescaleBranchLengths (0.5,1);
ExecuteCommands (constructLF ("codonLF", "codon_filter", "codon_tree", fileCount));

grid = defineAlphaBetaGrid (grid_points);

fprintf         (stdout, "[FUBAR PHASE 2] Determining appropriate branch scaling using the ", grid_points, "X", grid_points, " grid points.\n");
bestPair        = computeLFOnGrid ("codonLF", grid, 0);
bestPair        = bestPair%2;
bestAlpha = bestPair [Rows(bestPair)-1][0];
bestBeta  = bestPair [Rows(bestPair)-1][1];

if (bestAlpha < 0.01) {
    bestAlpha = 0.01;
}


LF_NEXUS_EXPORT_EXTRA = "";
rescaleBranchLengths (bestBeta/bestAlpha,0);
Export(lfExport,codonLF);
fprintf (codonFit, CLEAR_FILE,lfExport); 


fprintf         (stdout, "\tBest scaling achieved for dN/dS = ", Format(bestBeta/bestAlpha,5,2), ".\n\tComputing site-by-site likelihoods at ", 
                                     grid_points, "X", grid_points, " grid points\n");

gridInfo        = computeLFOnGrid ("codonLF", grid, 1);

fprintf         (stdout, "\tFinished with likelihood calculations. Achieved throughput of  ",
                                   Format(Rows(grid)/(Time(1)-timer),4,2), " calculations/second\n");


fprintf         (gridInfoFile,CLEAR_FILE, grid, "\n", gridInfo);

//------------------------------------------------------------------------------------------------//

function rescaleBranchLengths (dNdS, firstPass) {
    nonSynRate:=dNdS*synRate;
    ExecuteCommands   ("FindRoot(cf,`mgBrLen`-3,synRate,0,1e10);");
    nonSynRate=synRate;
    
    scalingFactor = cf/nf;
    
    global alpha = 1;
    global beta  = dNdS;
        
    for (file_part = 1; file_part <= fileCount; file_part += 1) {
        if (firstPass) {
            treeString = Eval ("Format (" + (nucLF_Info["Trees"])[file_part-1] + ",1,1)");
            ExecuteCommands   ("Tree codon_tree_" + file_part + " = " + treeString);
            ExecuteCommands   ("DataSetFilter codon_filter_" + file_part + " = CreateFilter (" + (nucLF_Info["Datafilters"])[file_part-1] + ",3,,,GeneticCodeExclusions)");
        } else {
            ExecuteCommands ("ClearConstraints (codon_tree_" + file_part + ");");
        }
        ExecuteCommands   ("ReplicateConstraint (\"this1.?.synRate:=alpha*scalingFactor__*this2.?.`paramName`__\",codon_tree_" + file_part + "," +  (nucLF_Info["Trees"])[file_part-1] + ");");
        ExecuteCommands   ("ReplicateConstraint (\"this1.?.nonSynRate:=beta*scalingFactor__*this2.?.`paramName`__\",codon_tree_" + file_part + "," +  (nucLF_Info["Trees"])[file_part-1] + ");");
    }
    return 0;
}

//------------------------------------------------------------------------------------------------//

function defineAlphaBetaGrid (one_d_points) {
    alphaBetaGrid = {one_d_points^2,2}; // (alpha, beta) pair
    oneDGrid      = {one_d_points,1};
   
    one_d_points    = Max (one_d_points, 10);
    neg_sel         = 0.7;
    neg_sel_points  = ((one_d_points)*neg_sel+0.5)$1;
    pos_sel_points  = (one_d_points-1)*(1-neg_sel)$1;
    if (neg_sel_points + pos_sel_points != one_d_points) {
        pos_sel_points = one_d_points - neg_sel_points; 
    }
    _neg_step = 1/neg_sel_points;
    for (_k = 0; _k < neg_sel_points; _k += 1) {
        oneDGrid [_k][0] =  _neg_step * _k;
    }
    oneDGrid [neg_sel_points-1][0] = 1;
    _pos_step = 49^(1/3)/pos_sel_points;
    for (_k = 1; _k <= pos_sel_points; _k += 1) {
        oneDGrid [neg_sel_points+_k-1][0] = 1+(_pos_step*_k)^3;
    }
    
    _p = 0;
    for (_r = 0; _r < one_d_points; _r += 1) {
        for (_c = 0; _c < one_d_points; _c += 1) {
           alphaBetaGrid[_p][0] = oneDGrid[_r];
           alphaBetaGrid[_p][1] = oneDGrid[_c];
           _p += 1;
        }
    }
    
    return alphaBetaGrid;   
}

fscanf              (stdin,"String", _sampleFile);
fscanf              (stdin,"String", _gridInfo);


fscanf              (stdin,"Number", _chainsToRun);
fscanf              (stdin,"Number", _chainLength);
fscanf              (stdin,"Number", _chainBurnin);
fscanf              (stdin,"Number", _chainSamples);
fscanf              (stdin,"Number", _concentration);

ExecuteAFile        (PATH_TO_CURRENT_BF + "FUBAR_tools.ibf");

assert (_chainsToRun > 1, "Must specify at least MCMC TWO chains to run");

/* the MCMC function */

baseFilePath  		= PATH_TO_CURRENT_BF + "spool/"+_in_FilePath;

debug = 0;

if (!debug) {
    funcText = "";
    funcsToExport = {"0": "runMCMC", "1": "jointLogL", "2": "LogDrichletDensity", "3": "computeLogLFfromGridAndWeights", "4": "siteLikelihoodsGivenWeights"};
    for (k = 0; k < Abs (funcsToExport); k+=1) {
        funcText += exportFunctionDefinition (funcsToExport[k]);
    }
    
    
     funcText += "\nfor (_chainIndex=start; _chainIndex<end; _chainIndex+=1) {runMCMC(_chainIndex,_gridInfo,_sampleFile,_chainLength,_chainBurnin,_chainSamples,_concentration);} return 0;";
    
     variablesToExport = "_gridInfo = \"" + _gridInfo + "\";\n" + 
                         "_sampleFile = \"" + _sampleFile + "\";\n" +   
                         "_chainsToRun = " + _chainsToRun + ";\n" +   
                         "_chainLength = " + _chainLength + ";\n" +   
                         "_chainBurnin = " + _chainBurnin + ";\n" +   
                         "_chainSamples = " + _chainSamples + ";\n" +   
                         "_concentration = " + _concentration + ";\n";

     if (MPI_NODE_COUNT > 1 && _chainsToRun > 1) {
            per_node    = Max(1,_chainsToRun $ MPI_NODE_COUNT);
            _startPoint = _chainsToRun-per_node;
            leftover    = _chainsToRun-per_node*MPI_NODE_COUNT;
            
            from          = 0;
            to            = per_node + (leftover>0);
            node_ranges   = {MPI_NODE_COUNT,2};
            
            for (node_id = 1; node_id < Min(_chainsToRun,MPI_NODE_COUNT); node_id += 1) {
                                        
                MPISend				(node_id, variablesToExport + ";start = " +from + ";end=" + to+";" + funcText); 
                
                
                node_ranges [node_id][0]         = from;
                node_ranges [node_id][1]         = to;
                
                from                             = to;
                to                              += per_node+(node_id<=leftover);  
            } 
        } else {
        _startPoint = 0;    
    }
        
    for (_r = _startPoint; _r < _chainsToRun; _r += 1){
        runMCMC(_r,_gridInfo,_sampleFile,_chainLength,_chainBurnin,_chainSamples,_concentration);
    }
    
    fprintf         (stdout, "\n[FUBAR PHASE 3 DONE] Finished running the MCMC chains; drew ", _chainsToRun, "x", _chainSamples, " samples from chains of length ", _chainLength, 
                             " after discarding ", _chainBurnin, " burn-in steps. Achieved throughput of ", Format(_chainLength/(Time(1)-time0),6,0) + " moves/sec.\n");
   
    if (MPI_NODE_COUNT > 1 && points > MPI_NODE_COUNT) {
        for (node_id = 1; node_id < Min(_chainsToRun,MPI_NODE_COUNT); node_id += 1) {
            MPIReceive (-1,fromNode,res);
        }
    }

    fprintf (_sampleFile,CLEAR_FILE, _chainsToRun, "\n");
}

//------------------------------------------------------------------------------------------------//

function jointLogL (weights, alpha) {
    ll  = computeLogLFfromGridAndWeights (weights);
    dir = LogDrichletDensity (weights, alpha);
    return {{ll__, dir__}};
}


//------------------------------------------------------------------------------------------------//
 

function LogDrichletDensity (dir_weights, alpha) {
     if (Min(dir_weights, 0) <= 1e-10) {
        return -1e10;
     }
     if (alpha == 1) {
        return 0;
     }  
     dim = Columns (dir_weights);
     return  (+dir_weights["Log(_MATRIX_ELEMENT_VALUE_)*(alpha-1)"]+LnGamma(alpha*dim)-dim*LnGamma(alpha));
}

//------------------------------------------------------------------------------------------------//

function computeLogLFfromGridAndWeights (wts) {
    return +(((wts *(gridInfo["conditionals"]))["Log(_MATRIX_ELEMENT_VALUE_)"])+gridInfo["scalers"]);
}

//------------------------------------------------------------------------------------------------//

function siteLikelihoodsGivenWeights (wts) {
    return wts*(gridInfo["conditionals"]);
}

//------------------------------------------------------------------------------------------------//

function runMCMC (chainID,gridFile, sampleFile, total,discard,expected_samples,_concentration_parameter){
    
    fscanf (gridFile, REWIND, "NMatrix,Raw", grid, gridInfo);
    gridInfo = Eval(gridInfo);
    
    points             = Rows(grid);
    sites              = Columns(gridInfo["conditionals"]);
    normalize_by_site  = ({1,points}["1"])*(gridInfo["conditionals"]);
    normalized_weights = (gridInfo["conditionals"])*({sites,sites}["1/normalize_by_site[_MATRIX_ELEMENT_ROW_]*(_MATRIX_ELEMENT_ROW_==_MATRIX_ELEMENT_COLUMN_)"]);
    sum_by_site        = normalized_weights * ({sites,1}["1"]);
   
    
    weights = {1,points}["Random(sum_by_site[_MATRIX_ELEMENT_COLUMN_]*0.8,sum_by_site[_MATRIX_ELEMENT_COLUMN_]*1.2)"];
    weights = weights * (1/(+weights));
    
    gridSampled = {points, 3};
    for (k = 0; k < points; k+=1) {
        gridSampled[k][0] = grid[k][0];
        gridSampled[k][1] = grid[k][1];
        gridSampled[k][2] = weights[k];
    }
    
    //fprintf (stdout, +weights["_MATRIX_ELEMENT_VALUE_<1e-10"], "\n");
    
    defaultStep = Max(Min(0.001,1/sites),(weights%0)[points*50$100]);
    
    //fprintf (stdout, "\nDefault step = ", defaultStep, "\n");
     
    currentSiteLikelihoods = siteLikelihoodsGivenWeights (weights);
    currentSiteLogSum      = +(currentSiteLikelihoods["Log(_MATRIX_ELEMENT_VALUE_)"]);
    currentLogL            = jointLogL (weights, _concentration_parameter);
    individualGridPointContribs = {};
    
    for (k = 0; k < points; k+=1) {
        individualGridPointContribs [k] = (gridInfo["conditionals"])[k][-1];
    }    
    
    contracting            = total*50;
    sample                 = (total-discard)$expected_samples;
    sampled_weights        = {expected_samples,points};
    sampled_likelihoods    = {1,expected_samples};
    
    time0                = Time(1);
    sample_index         = 0;
    
    baselineStep         = defaultStep;
    reductionFactor      = 1;
    accepted_steps       = 0;

    if (MPI_NODE_ID == 0) {
        fprintf         (stdout, "\n[FUBAR PHASE 3] Running an MCMC chain (ID ",chainID, ") to obtain a posterior sample of grid point weights: ", total, 
                                            " total steps, of which ", discard, " will be discarded as burn-in, and sampling every ", sample, " steps. Dirichlet prior concentration parameter = ", _concentration, ".\n");
        if (MPI_NODE_COUNT > 1)
        {
            fprintf         (stdout, "\tIn addition, ", _chainsToRun-1 , " independent chains (started at random points in the search spaces) are being run in parallel to evaluate convergence and effective sample sizes. The final sample will include thinned post burn-in representatives from all chains\n");
        }
    }

    totalStepSum        = 0;
    meanSampledLogL     = 0;
    initialLogL         = currentLogL[0];

    for (steps = 0; steps < total; steps += 1) {              

        idx    = Random (0, points-1e-10)$1;
        idx2   = Random (0, points-1e-10)$1;
        while (idx == idx2) {
            idx2   = Random (0, points-1e-10)$1;     
        }
        
        if ((steps+1) % contracting == 0) {
            acc_rate = accepted_steps/steps;
            if (acc_rate < 0.25) {
                baselineStep = baselineStep/1.6;
            } else if (acc_rate > 0.5) {
                baselineStep = baselineStep*1.6;                
            }
        }
        
        change = Random (0,baselineStep);
        totalStepSum += change;
        
        if (weights[idx] > change) {
            diffVector          = (individualGridPointContribs[idx2]-individualGridPointContribs[idx])*change;
            logLDiff            = +((currentSiteLikelihoods + diffVector)["Log(_MATRIX_ELEMENT_VALUE_)"]) - currentSiteLogSum;
            diffPrior           = (_concentration_parameter-1)*(Log((weights[idx]-change)/weights[idx])+Log((weights[idx2]+change)/weights[idx2]));
            costOfMove          = logLDiff+diffPrior;
            
            if (Random (0,1) <= Exp (costOfMove)) {
            
                currentLogL[0] += logLDiff;
                currentLogL[1] += diffPrior;
         
                currentSiteLikelihoods += diffVector;
                currentSiteLogSum += logLDiff;
                 
                weights[idx]  += (-change);
                weights[idx2] += (+change);
                accepted_steps += 1;
            } 
        }
        
        if (steps > discard) {
            if ((steps - discard + 1) % sample == 0) {
                for (dd = 0; dd < points; dd += 1) {
                    sampled_weights[sample_index][dd] = weights[dd];
                }
                sampled_likelihoods[sample_index] = currentLogL[0];
                meanSampledLogL += currentLogL[0];
                sample_index += 1;
                logLString = meanSampledLogL/sample_index;
            }
        } else {
            logLString = "(burning in)";
        }
    
        if ((1+steps) % sample == 0) {
             if (MPI_NODE_ID == 0) {
                SetParameter (STATUS_BAR_STATUS_STRING, "Running MCMC chain ID "+ chainID + ". Current step: " + (1+steps) + "/" + total + ". Mean sampled log(L) = " + logLString 
                + ". Acceptance rate = " + accepted_steps/steps, 
                0);
            }
        }
    }

    mcmcfile = _sampleFile + "." + chainID;
    fprintf (mcmcfile,CLEAR_FILE, sampled_likelihoods, "\n\n", sampled_weights);
    
    return 0;
}
fscanf              (stdin, "String", nuc_fit_file);
fscanf              (stdin, "String", grid_file);
fscanf              (stdin, "String", sample_base_file);
fscanf              (stdin, "Number", _chainCount);
fscanf              (stdin, "String", results_file);

// for PHASE 5
fscanf (stdin, "String", sim_fit_file);
fscanf (stdin, "String", sim_grid_info);
fscanf (stdin, "String", codon_fit_file);

ExecuteAFile        (PATH_TO_CURRENT_BF + "FUBAR_tools.ibf");
LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("WriteDelimitedFiles");

ExecuteAFile        (nuc_fit_file);
sequences = nucData_1.species;
GetInformation (treeCount,"^nuc_tree_[0-9]+$");
fileCount       = Columns (treeCount);
treeLengths     = {fileCount,1};

for (fileID = 1; fileID <= fileCount; fileID += 1)
{
	treeLengths [fileID-1] = + Eval("BranchLength(nuc_tree_"+fileID+",-1)");
}

fscanf (grid_file, REWIND, "NMatrix,Raw", grid, site_probs);
site_probs = Eval (site_probs);
sites   = Columns (site_probs["conditionals"]);


readMCMCSamples (sample_base_file,_chainCount);


notPositiveSelection = {points,1} ["grid[_MATRIX_ELEMENT_ROW_][0]>=grid[_MATRIX_ELEMENT_ROW_][1]"];
nonPositiveCount     = +notPositiveSelection;

priorMean            = {1, points};
sampleFromThisDistro = {nonPositiveCount,2};

tabulateGridResults (points, sites, samples, _chainCount);

from = 0;
for (_point = 0; _point < points; _point += 1) {
    priorMean [_point] = (+jointSamples[-1][_point])/samples;
    if (notPositiveSelection [_point]) {
        sampleFromThisDistro [from][0] = _point;
        sampleFromThisDistro [from][1] = priorMean [_point];
        from += 1;
    }
}

priorNN = +(sampleFromThisDistro [-1][1]);
fubar_results = reportSiteResults   (sites, 0, priorNN, _fubar_do_simulations);
fubarRowCount     = Rows (fubar_results);

if (_fubar_do_simulations) {
    simPatterns = Random(sampleFromThisDistro, {"PDF":"Multinomial","ARG0":1000});
    fprintf (sample_base_file, CLEAR_FILE, _chainCount, "\n", jointLogL, "\n", jointSamples);
    
    fprintf (stdout, "\n[FUBAR PHASE 5] Performing 1,000 simulations under the data-derived composite null model to derive False Discovery Rate (FDR) estimates\n");
    ExecuteAFile (PATH_TO_CURRENT_BF + "FUBAR_PHASE_5.bf");
    
    // sort posteriorsUnderNN on the posterior prob of pos.sel. column
    
    posteriorsNNPP    = (posteriorsUnderNN[-1][3]) % 0;
    sortedFubarP      = ({Rows (fubar_results), 2} ["_MATRIX_ELEMENT_ROW_*(_MATRIX_ELEMENT_COLUMN_==0)+fubar_results[_MATRIX_ELEMENT_ROW_][3]*(_MATRIX_ELEMENT_COLUMN_==1)"])%1;
    
    currentNNRows     = Rows (posteriorsUnderNN);
    
    currentNNIndex    = currentNNRows - 1;
    
    
    for (currentFubarIndex =  fubarRowCount - 1; currentFubarIndex >= 0; currentFubarIndex += -1) {
        currentFubarPosteriorP = sortedFubarP[currentFubarIndex][1];
        
        while (currentNNIndex > 0 && posteriorsNNPP[currentNNIndex] > currentFubarPosteriorP) {
            currentNNIndex = currentNNIndex - 1;
        }
        
        FDR = Min (1,(currentNNRows-currentNNIndex)/currentNNRows * priorNN / ((fubarRowCount-currentFubarIndex) / fubarRowCount));
        
        if (currentNNIndex == 0) {
            break;
        }
        fubar_results[sortedFubarP[currentFubarIndex][0]][7] = FDR;
        
    }
    
    for (; currentFubarIndex >= 0; currentFubarIndex += -1) {
        FDR = Min (1, priorNN / ((fubarRowCount-currentFubarIndex) / fubarRowCount));
        fubar_results[sortedFubarP[currentFubarIndex][0]][7] = FDR;
    }
}

site_counter = {};
for (currentFubarIndex = 0; currentFubarIndex < fubarRowCount; currentFubarIndex += 1) {
    site_counter + (currentFubarIndex+1);
}

if (_fubar_do_simulations) {
    WriteSeparatedTable (results_file, {{"Codon","alpha","beta","beta-alpha","Prob[alpha<beta]", "Prob[alpha>beta]", "BayesFactor","PSRF", "Neff", "FDR"}}, fubar_results, site_counter, ",");
} else {
    WriteSeparatedTable (results_file, {{"Codon","alpha","beta","beta-alpha","Prob[alpha<beta]", "Prob[alpha>beta]", "BayesFactor","PSRF", "Neff"}}, fubar_results, site_counter, ",");
}
 ExecuteAFile (codon_fit_file);

codonCharactersAgrument = {{"A","C","G","T"}
			  			   {"3",GeneticCodeExclusions,"",""}};

nn_sample_count = + (simPatterns[-1][1]);
byGridPoint     = {nn_sample_count, 3};

GetString         (lfInfo, codonLF, -1);
_treeCount      = Columns (lfInfo["Trees"]);
_filterNames    = lfInfo ["Datafilters"];

baseFrequencies = Eval ((lfInfo["Base frequencies"])[0]);
GetString (targetSpeciesOrdering, existingFit.ds_1,-1);

sequenceCount = existingFit.ds_1.species;
indexer       = 0;


accumulatedSequences = {};

for (k = 0; k < _treeCount; k+=1) {
    accumulatedSequences [k] = {};
    for (k2 = 0; k2 < sequenceCount; k2+=1) {
        (accumulatedSequences [k])[k2] = "";
        (accumulatedSequences [k])[k2] *128;
    }
}

for (k = 0; k < Rows(simPatterns); k += 1) {
    samples_for_this_value = simPatterns[k][1];
    alpha_beta  = grid[simPatterns[k][0]][-1];
    for (k2 = 0; k2 < samples_for_this_value; k2+=1) {
        byGridPoint[indexer][0] = alpha_beta[0];
        byGridPoint[indexer][1] = alpha_beta[1];
        byGridPoint[indexer][2] = 1+Random(0,_treeCount-0.000000001)$1;
        indexer += 1;
    }
}

byGridPoint  = byGridPoint % {{0,1,2}};

//fprintf (sim_grid_info, CLEAR_FILE, byGridPoint);

extraFunctions        = exportFunctionDefinition("simulateFromTreeGivenRates") + exportFunctionDefinition("mapSets");
LF_NEXUS_EXPORT_EXTRA = "\n;PRESERVE_SLAVE_NODE_STATE=1;MPI_NEXUS_FILE_RETURN=Rows(UserFunction);\nbaseFrequencies= " + baseFrequencies + ";\ncodonCharactersAgrument="+codonCharactersAgrument+"\n;targetSpeciesOrdering="+targetSpeciesOrdering+";\n"; 
                        

MPI_NODE_STATUS = {MPI_NODE_COUNT,1};

Export(lfExport,codonLF);

for (k = 1; k < MPI_NODE_COUNT; k+=1) {
    MPISend (k, lfExport);
}

for (k = 1; k < MPI_NODE_COUNT; k+=1) {
    MPIReceive (-1, fromNode, result);
 }

//------------------------------------------------------------------------------

currentIndex = 0;
t0 = Time(1);
for (k = 1; k < nn_sample_count; k+=1)
{
    if (Abs(byGridPoint [k][-1] - byGridPoint[k-1][-1]) > 1e-10) {
        treeID         = byGridPoint[k-1][2];
        dispatchABlock (byGridPoint[k-1][2],byGridPoint[k-1][0], byGridPoint[k-1][1],  k-currentIndex);
        currentIndex = k;
        SetParameter (STATUS_BAR_STATUS_STRING, "Simulating neg/neutral sites "+ (k) + "/" + nn_sample_count + " " + _formatTimeString(Time(1)-t0),0);
    }
}
   
dispatchABlock (byGridPoint[k-1][2],byGridPoint[k-1][0], byGridPoint[k-1][1],  k-currentIndex);

still_busy = +(MPI_NODE_STATUS["_MATRIX_ELEMENT_VALUE_>0"]);
for (k = 0; k < still_busy; k += 1) {
    processABlock (0);
}


for (k = 0; k < _treeCount; k+=1) {
    newData = ""; newData * 128;
    for (k2 = 0; k2 < sequenceCount; k2+=1) {
         (accumulatedSequences [k])[k2] *0;
         newData * (">" + targetSpeciesOrdering[k2] + "\n" + (accumulatedSequences [k])[k2] + "\n");
    }
    newData * 0;
    ExecuteCommands ("DataSet simulated_data_" + (k+1)  +  " = ReadFromString (newData)");
    ExecuteCommands ("GetDataInfo (info, " + _filterNames[k] +",\"PARAMETERS\");");
    ExecuteCommands ("DataSetFilter " + _filterNames[k] + " = CreateFilter (simulated_data_" + (k+1) + "," + info["ATOM_SIZE"] + ",,\"" + info["SEQUENCES_STRING"] + "\",\"" + info["EXCLUSIONS"] + "\");");
}

LIKELIHOOD_FUNCTION_OUTPUT = 7;
LF_NEXUS_EXPORT_EXTRA = "";
fprintf (sim_fit_file, CLEAR_FILE, codonLF);
DeleteObject (codonLF);
ExecuteAFile (sim_fit_file);
site_probs = computeLFOnGrid ("codonLF", grid,1);

fprintf (sim_grid_info, CLEAR_FILE, byGridPoint, "\n", site_probs);

tabulateGridResults (points, nn_sample_count, samples, _chainCount);
posteriorsUnderNN = reportSiteResults (nn_sample_count, 0, 0, 1);

return posteriorsUnderNN;

//------------------------------------------------------------------------------

function   dispatchABlock (treeID, alphaV, betaV, howMany) {
    if (MPI_NODE_COUNT <= 1) {
        simulatedChunk = simulateFromTreeGivenRates (treeID,alphaV,betaV,howMany, "baseFrequencies", "codonCharactersAgrument", targetSpeciesOrdering);
        return processABlock (treeID);
    }

    //fprintf(stdout, MPI_NODE_STATUS, "\n");
    
    for (node_id = 1; node_id < MPI_NODE_COUNT; node_id += 1) {
        if (MPI_NODE_STATUS[node_id] == 0) {
            break;
        }
    }
    //fprintf (stdout, "Node id:", node_id,"\n");
    
    if (node_id == MPI_NODE_COUNT) {
        node_id = processABlock(0);
    }
    
    mpi_command = extraFunctions + "\nreturn simulateFromTreeGivenRates ("+treeID+","+alphaV+","+betaV+","+howMany+",\"baseFrequencies\", \"codonCharactersAgrument\", targetSpeciesOrdering);";
    //fprintf (stdout, "[DEBUG:] To node ", node_id, "=>", mpi_command, "\n");
    MPISend (node_id, mpi_command);
    //fprintf (stdout, MPI_LAST_SENT_MSG);
    MPI_NODE_STATUS[node_id] = treeID;
    return 0;
}

//------------------------------------------------------------------------------

function   processABlock (treeID) {
   if (MPI_NODE_COUNT > 1) {
        MPIReceive (-1,fromNode,result);
        simulatedChunk = Eval (result);
        //fprintf (stdout, "Result:",result,"\n");
        treeID = MPI_NODE_STATUS[fromNode];
        MPI_NODE_STATUS[fromNode] = 0;
   }    
   
   for (k2 = 0; k2 < sequenceCount; k2+=1) {
        ((accumulatedSequences [treeID-1])[k2]) * simulatedChunk[k2];
   }           
   
   return fromNode;
}

//------------------------------------------------------------------------------

function simulateFromTreeGivenRates (treeID, alphaValue, betaValue, codons, frequencies, characters, speciesOrdering) {
    alpha = alphaValue;
    beta  = betaValue;
    ExecuteCommands ("DataSet theData = Simulate (codon_tree_"+treeID+",`frequencies`,`characters`,codons,0);");
    GetString (similatedOrdering, theData,-1);
    reordering = mapSets (speciesOrdering, similatedOrdering);
    DataSetFilter reporterFilter = CreateFilter (theData,1,"",Join(",",reordering));
    GetInformation (allSeqs, reporterFilter);
    GetString (similatedOrdering, reporterFilter,-1);
    return allSeqs;
}
