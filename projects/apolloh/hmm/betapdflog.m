function y = betapdflog(x,a,b)
%% FUNCTION betapdflog(x,alphaParam,betaParam)
%
% author: Gavin Ha <gha@bccrc.ca>
%         Dept of Molecular Oncolgy/Centre for Translational and Applied Genomics
%         British Columbia Cancer Agency
%         University of British Columbia
% date  : April 11, 2011

y = -betaln(a,b) + (a-1)*log(x) + (b-1)*log(1-x);