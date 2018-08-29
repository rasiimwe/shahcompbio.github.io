function [mu_N] = estimateBetaParamsMap(alpha,beta)
%% FUNCTION estimateBetaParamsMap(alpha,beta)
%
% author: Gavin Ha <gha@bccrc.ca>%         
%         Dept of Molecular Oncolgy/Centre for Translational and Applied Genomics
%         British Columbia Cancer Agency
%         University of British Columbia
% date  : April 10, 2011    
% given the previous estimate of the parameters and the current 
% estimate of the reponsibilities, estimate the new parameters
%

% alpha  - parameter to Beta prior, 1xK
% beta   - parameter to Beta prior, 1xK
%
% RETURN:
% =======
% mu_N     - new MAP estimate of means

 
% mu_N = (alpha - 1)/(alpha + beta - 2) 
mu_N = (alpha - 1)./(alpha + beta - 2);