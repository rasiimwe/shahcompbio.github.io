function rho = distanceTransitionFunction(d,L)
% FUNCTION rho = distanceTransitionFunction(d,L)
%
% author: Gavin Ha <gha@bccrc.ca>
%         Dept of Molecular Oncolgy/Centre for Translational and Applied Genomics
%         British Columbia Cancer Agency
%         University of British Columbia
% date  : March 7, 2011
d = double(d);
rho = (1/2)*(1-exp(-d/(2*L)));