function f = stromalDerivativeUpdateEqn(s,rt,rn,a,b,alphaS,betaS)
%% FUNCTION binomialpdf(k,N,mu)
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
% date  : August 26, 2011

K = length(rt);

% data likelihood derivative wrt to s
dlik_ds = 0;
for k=1:K
    dlik_ds = dlik_ds + (rn-rt(k)) * ((a(k)/(s*rn+(1-s)*rt(k))) - (b(k)/(1-s*rn-(1-s)*rt(k))));
end

% beta prior likelihood (of s) derivative wrt to s
dbeta_ds = (alphaS-1)/s - (betaS-1)/(1-s);

f = dlik_ds + dbeta_ds;
