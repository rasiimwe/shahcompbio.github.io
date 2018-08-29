function [s_new, rt_new, pi] = estimateBinomNoiseStromalRatioParamsMap(x,N,rho,estimateS,rt,rn,alphaR,betaR,alphaS,betaS,kappa)
%% FUNCTION estimateBinomNoiseParamsMap(x,N,rho,alpha,beta,kappa
%
% author: Gavin Ha <gha@bccrc.ca>
%         Sohrab Shah <sshah@cs.ubc.ca>
%         Dept of Molecular Oncolgy/Centre for Translational and Applied Genomics
%         British Columbia Cancer Agency
%         University of British Columbia
% date  : March 7, 2011    
% given the previous estimate of the parameters and the current 
% estimate of the reponsibilities, estimate the new parameters
%
% x      - data, 1xT
% N      - Number of trials, 1xT
% rho    - current estimate of the responsibilites, KxT
% alpha  - parameter to Beta prior, 1xK
% beta   - parameter to Beta prior, 1xK
% kappa  - prior to stationary distribution, 1xK
% RETURN:
% =======
% mu_N     - new MAP estimate of means
% pi       - new MAP estimate of stationary distribution

    
[K,T] = size(rho);
rt_new = zeros(1,K);
% vectorize all parameters

% stromal parameter, Beta prior
a = x * rho';
b = (N-x) * rho';
interval = [eps,1-eps];
if estimateS==1  %if we want to estimate S
    s_new = fzero(@(s) stromalDerivativeUpdateEqn(s,rt,rn,a,b,alphaS,betaS),interval);
else
    s_new = 0;
end

% tumour reference allelic ratio parameter, Beta prior

for k=1:K
    rt_new(k) = fzero(@(r) allelicRatioDerivativeUpdateEqn(s_new,r,rn,a(k),b(k),alphaR(k),betaR(k)),interval);
end
rt_new = rt_new';

% stationary distribution, Dirichlet prior
% pi(k) = (kappa(k) + sum(rho(k)) - 1)/(sum(kappa) + N - K)
pi = (rho(:,1)+kappa'-1)./(nansum(rho(:,1)) + nansum(kappa) - K);
pi = pi';