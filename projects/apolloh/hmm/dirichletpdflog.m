function y = dirichletpdflog(x,k)
%% FUNCTION dirichletpdflog(x,k)
%
% author: Gavin Ha <gha@bccrc.ca>
%         Dept of Molecular Oncolgy/Centre for Translational and Applied Genomics
%         British Columbia Cancer Agency
%         University of British Columbia
% date  : April 11, 2011

%normalizing constant
c = gammaln(sum(k)) - sum(gammaln(k));

%likelihood
l = sum((k-1).*log(x));

%together
y = c + l;