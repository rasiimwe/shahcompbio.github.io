function [pi] = estimateDirichletParamsMap(kappa)
%% FUNCTION estimateBinomNoiseParamsMap(x,N,rho,alpha,beta,kappa
%
% author: Gavin Ha <gha@bccrc.ca>
%         Dept of Molecular Oncolgy/Centre for Translational and Applied Genomics
%         British Columbia Cancer Agency
%         University of British Columbia
% date  : March 7, 2011    
% given the previous estimate of the parameters and the current 
% estimate of the reponsibilities, estimate the new parameters
%
% kappa  - prior to stationary distribution, 1xK
%
% RETURN:
% =======
% pi       - new MAP estimate of stationary distribution

[jnk,K] = size(kappa);

% calculate the stationary distribution
% pi(k) = (kappa(k) + sum(rho(k)) - 1)/(sum(kappa) + N - K)
pi = (kappa-1)./(nansum(kappa) - K);