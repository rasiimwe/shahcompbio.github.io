function Z = apolloh(infile,cnfile,paramset,outparam,outfile)
% Function apolloh(infile,cnfile,paramset,outparam,outfile)
% 
% APOLLOH v0.1.0
% Hidden Markov model (HMM) for predicting somatic loss of heterozygosity
% and allelic imbalance in tumour whole genome sequencing data
%
% author: Gavin Ha <gha@bccrc.ca>
%         Dept of Molecular Oncolgy/Centre for Translational and Applied Genomics
%         British Columbia Cancer Agency
%         University of British Columbia
% date  : December 15, 2011
%
% INPUTS:
%
% infile         Tab-delimited input file containing allelic counts from
%                the tumour at positions determined as heterozygous from
%                the normal genome.
%                6 columns: 1) chr 2) position 3) reference base 4)
%                referenc count 5) non-reference base 6) non-reference
%                count 
%
% cnfile         Tab-delimited input copy number segment prior file.
%                The accepted format is the output from HMMcopy, 
%                a read-depth for analyzing copy number in tumour-
%                normal sequenced genomes.  
%                However, copy number segments from any source can be used.
%                8-columns:
%                 1) id (can be arbitrary, not used)
%		  2) chr
%                 3) start  
%		  4) stop
%                 5) Number of 1kb intervals (can be arbitrary; not used)
%		  6) median log2 ratio (normal and tumour) for segment
%                 7) HMM state: 1=HOMD (0 copies), 2=HEMD (1 copy),
%                        3=NEUT (2 copies), 4=GAIN (3 copies),
%                        5=AMP (4 copies), 6=HLAMP (5+ copies)
%		         Note that for AMP and HLAMP, relative numbers of copies
%                        can be used (i.e. GAIN is 3-4 copies, AMP is 5-6 copies,
%	         	HLAMP is 7+ copies)
%		  8 )CN state (can be arbitrary; not used)
%			If cnfile='0' is used, then copy number of 2 (diploid) is used
%                       for all positions.   
%
% paramset       Parameter intialization file is a matlab binary (.mat) file.  
%                This file contains model and setting paramters necessary 
%                to run the program.  
%
%                See examples in "<install_dir>/parameters/".
%
% OUTPUTS:
%
% outfile        Tab-delimited output file for position-level results. 
%                  9-columns:
%                     1) chr ('X' and 'Y' will be output as 23 and 24)
%                     2) position
%                     3) reference count
%                     4) non-reference count
%                     5) total depth
%                     6) allelic ratio
%                     7) copy number (from input)
%                     8 ) APOLLOH genotype state
%                     9) Zygosity state.
%                N additional columns: posterior marginal probabilities
%                (responsibilities) for each APOLLOH genotype state.
%                Zygosity states are: 
%                DLOH=deletion-LOH (state 1) 
%                NLOH=copy-neutral-LOH (states 2,4) 
%                ALOH=amplified-LOH (states 5,8,9,13,14,19) 
%                HET=heterozygous (states 3,6,7)
%                ASCNA=allele-specific-amplification (states 10,12,15,18) 
%                BCNA=balanced-amplification (states 11,16,17)
%                
%                Segment boundaries are determined as consecutive marginal
%                states of DLOH, NLOH, ALOH, HET, BCNA, ASCNA; this 
%                implementation does not output this information. An 
%                external Perl script handles this: 
%                "<install_dir>/scripts/analysis/createSingleSegFileFromAPOLLOH.pl" 
%
% outparam       Tab-delimited output file storing converged parameters after 
%                model training using Expectation Maximization (EM) algorithm.  
%                1) Number of iterations
%                2) Global normal contamination parameter
%                3) Binomial parameters for each HMM class/state. 
%

disp('Running APOLLOH'); 
%% Read in data
disp(['apolloh: Loading data ',infile]);
% Genome data
%1  799550        G  48    C     3     0.9999999983    0.0000000017    0.0000000000    1       aa
fid = fopen(infile);
data = textscan(fid,'%s%d%s%d%s%d','delimiter','\t');
chr = data{1}; chr(strcmp(chr,'X')) = {'23'}; chr(strcmp(chr,'Y')) = {'24'}; chr(strcmp(chr,'MT')) = {'25'};
chr = str2double(chr);
posn = double(data{2});
ref = double(data{4});
nonRef = double(data{6});
depth = ref + nonRef;

% load parameters from a script
disp(['apolloh: Loading parameters from ',paramset]);

%% Load parameters from configuration file
%   contains the following parameters:
%   alphaHyperS, betaHyperS, alphaHyperR, betaHyperR, kappaHyper, maxDepth, minDepth, mirror, txnLen, maxiter,
%   rt,rn,
%
load(paramset);

%% Filter
keepNormalDepth = depth <= maxDepth & depth >= minDepth;

%% get copy number for each position in data; if no cnfile given, then use copy number of 2 for all data points
if ~strcmp(cnfile,'0')
    cn = getCNfromSegs(chr,posn,cnfile,11,0);
else
    cn = repmat(2,[length(ref),1]);
end
nonHOMD = cn > 0;
cI = chr<=23 & keepNormalDepth & nonHOMD;

%% Set up and run HMM
T = length(ref(cI));
if (mirror == 1)
    refToUse = max(ref,nonRef);
    CNS = [0,1,3,5,8,11];
    ZS = [1,2,1,2,1,3,2,1,3,2];
else
    refToUse = ref;
    CNS = [0,1,4,8,13,19];
    ZS = [1,2,1,1,2,2,1,1,3,2,3,1,1,3,2,2,3,1];
end
ratio = refToUse ./ depth;

% initialize model parameters
if estimateS==1
    s_0 = estimateBetaParamsMap(alphaHyperS,betaHyperS);
else
    s_0 = 0;
end
r_0 = estimateBetaParamsMap(alphaHyperR,betaHyperR);
pi_0 = estimateDirichletParamsMap(kappaHyper);

%%%%%%%%%%%%%%%%%%%% run Forward-Backward and Viterbi %%%%%%%%%%%%%%%%%%%
[partG,s,rt,mus,pi,loglik,partRho,segs] = binomHMMfwdBackViterbiCStromalRatio(refToUse(cI),depth(cI),chr(cI),posn(cI),cn(cI),r_0,rn,estimateS,s_0,alphaHyperR,betaHyperR,alphaHyperS,betaHyperS,pi_0,kappaHyper,txnExpLen,ZS,CNS,maxiter);

G = zeros(length(ref),1)-1; %full list of points
G(cI) = partG; %assign analyzed points
G(~nonHOMD) = 0; %use state 0 for homd
rho = zeros(length(r_0),length(ref));
rho(:,cI) = partRho; %homd positions will have 0 posterior probs
clear partRho;
%% OUTPUT
%%%%%%%%%%%%%%%%%%%% output Z (for each datapoint) %%%%%%%%%%%%%%%%%%%
outInd = chr<=23 & keepNormalDepth; %include HOMD in output
%re-annotate with deletion cn and DLOH
if ~strcmp(cnfile,'0')
    cnWithDel = getCNfromSegs(chr,posn,cnfile,11,1); 
else
    cnWithDel = repmat(2,[length(ref),1]);
end
G(outInd & nonHOMD) = G(outInd & nonHOMD) + 1; %add 1 to include DLOH as state 1
if mirror==1 
    dlohInd = G==2 & cnWithDel==1;  %check for state LOH (neutral) and copy number del of 1
else
    dlohInd = (G==2 | G==4) & cnWithDel==1;
end
G(dlohInd) = 1;  %assign DLOH as state 1
Z = decodeLOH(G(outInd),mirror);
rho = rho(:,outInd);
disp(['apolloh: Writing results to ',outfile]);
fid = fopen(outfile,'w+');
formatStr = '';
fullmatrix = [chr(outInd),posn(outInd),ref(outInd),depth(outInd)-ref(outInd),depth(outInd),ref(outInd)./depth(outInd),cnWithDel(outInd),G(outInd)];
for k=1:length(r_0)
    formatStr = [formatStr,'\t%1.4f'];
end
for i=1:length(Z)
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%1.4f\t%d\t%d',fullmatrix(i,:));
    fprintf(fid,'\t%s',Z{i});        
    fprintf(fid,[formatStr,'\n'],rho(:,i)');   
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%% output segments of Z %%%%%%%%%%%%%%%%%%%
%disp(['apolloh: Writing segment results to ',segoutfile]);
%write out all segments
%outputSegsToFile(segs,chr(cI),posn(cI),ratio(cI),mirror,segoutfile)
%%%%%%%%%%%%%%%%%%%% save data to mat file %%%%%%%%%%%%%%%%%%%
disp(['apolloh: Saving results to ',outparam]);
fid = fopen(outparam,'w+');
fprintf(fid,'Number of iterations:\t%d\n',int32(length(s)));
fprintf(fid,'Normal cell contamination:\t%0.4f\n',s(end));
mus_str = sprintf('%0.4f ' ,mus(:,end));
mus_str = strtrim(mus_str);
fprintf(fid,'Binomial Means:\t%s\n',mus_str);
fclose(fid);
%save(paramfile,'G','Z','mus','pi','loglik','s','rt');
%%%%%%%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%%%%%%
%disp('apolloh: Plotting results...');
%plotResponsibilities(rho,chr(cI),posn(cI));
%plotLogLikelihood(loglik);
%plotSimulatedDepthPoisson('',depth(cI),cn(cI));
%plotParamEM(mus);
% [path,name,ext] = fileparts(outfile);
% mkdir([path,name]);
% outplot = [path,'/',name,'/',name];
% if (mirror==1)
%     plotBinoPDFZ(0:40,40,mus(:,end),outplot);
% else
%     plotBinoPDF(0:40,40,mus(:,end),outplot);
% end
% plotSNPallelicRatio(outfile,'',outplot);
