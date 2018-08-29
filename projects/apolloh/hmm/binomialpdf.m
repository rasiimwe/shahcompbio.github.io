function y = binomialpdf(k,N,mu)
%% FUNCTION binomialpdf(k,N,mu)
%
%  Computed in log form but returned in exponentialed form
%
% author: Gavin Ha <gha@bccrc.ca>
%         Dept of Molecular Oncolgy/Centre for Translational and Applied Genomics
%         British Columbia Cancer Agency
%         University of British Columbia
% date  : April 11, 2011

%normalizing constant
c = gammaln(N+1) - gammaln(k+1) - gammaln(N-k+1);

%likelihood
l = k*log(mu) + (N-k)*log(1-mu);

%together
y = c + l;
y = exp(y);