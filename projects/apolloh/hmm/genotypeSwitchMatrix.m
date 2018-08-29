function C = genotypeSwitchMatrix(C,cn,CNS)
%% FUNCTION CY = genotypeSwitchMatrix(K,cn,CNS)
%
% author: Gavin Ha <gha@bccrc.ca>
%         Dept of Molecular Oncolgy/Centre for Translational and Applied Genomics
%         British Columbia Cancer Agency
%         University of British Columbia
% date  : April 12, 2011
%
% cn = copy number data for each position

[K,N] = size(C);

for i=1:N
   C(CNS(cn):CNS(cn+1),i) = 1; 
end
