function f = allelicRatioDerivativeUpdateEqn(s,rt,rn,a,b,alphaR,betaR)
%% FUNCTION binomialpdf(k,N,mu)
%
%  2-component mixture for the Binomial mean given stromal contamination
%  parameter (s) and theoretical means (r)
%  This function evaluates partial derivative for only one state
%
%  s = stromal contamination parameter
%  rt = scalar allelic ratio for state k; theoretical reference allelic ratio; 
%       values should reflect any skew in data 
%  rn = scalar representing reference allelic ratio of heterozygous normal
%       sample; should also reflect any skew in data, i.e. 0.5 for no skew,
%       0.6 for skew of data towards reference
%  alphaR, betaR = scalar hyperparameter for state k 
%
% author: Gavin Ha <gha@bccrc.ca>
%         Dept of Molecular Oncolgy/Centre for Translational and Applied Genomics
%         British Columbia Cancer Agency
%         University of British Columbia
% date  : August 26, 2011

% data likelihood derivative wrt to s
dlik_dr = (1-s) * ((a/(s*rn+(1-s)*rt)) - (b/(1-s*rn-(1-s)*rt)));

% beta prior likelihood (of s) derivative wrt to s
dbeta_dr = (alphaR-1)/rt - (betaR-1)/(1-rt);

f = dlik_dr + dbeta_dr;
