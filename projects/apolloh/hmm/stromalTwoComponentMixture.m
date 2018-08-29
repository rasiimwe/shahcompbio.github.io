function mus = stromalTwoComponentMixture(rt,rn,s)
%% FUNCTION mus = stromalTwoComponentMixture(rt,rn,s)
%
%  2-component mixture for the Binomial mean given stromal contamination
%  parameter (s) and theoretical means (r)
%
%  s = stromal contamination parameter
%  rt = vector of length K (# of states); theoretical reference allelic ratio; 
%       values should reflect any skew in data 
%  rn = float representing reference allelic ratio of heterozygous normal
%       sample; should also reflect any skew in data, i.e. 0.5 for no skew,
%       0.6 for skew of data towards reference
%
% author: Gavin Ha <gha@bccrc.ca>
%         Dept of Molecular Oncolgy/Centre for Translational and Applied Genomics
%         British Columbia Cancer Agency
%         University of British Columbia
% date  : August 23, 2011

mus = s*rn + (1-s)*rt;
