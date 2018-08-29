function [gamma,alpha,beta,loglik] = nonStationaryFwdBack(init_state_distrib,transmat,obslik)

% forwards propagation, backwards sampling
%
% input
% init_state_distrib(q)
% transmat(q, q')
% obslik(q, t)
%
% output


[Q T] = size(obslik);
scale = ones(1,T);

alpha = zeros(Q,T);
beta = zeros(Q,T);
gamma = zeros(Q,T);

t = 1;
alpha(:,1) = init_state_distrib(:) .* obslik(:,t);
[alpha(:,t), scale(t)] = normalise(alpha(:,t));
for t=2:T
    trans = reshape(transmat(:,:,t),[Q,Q]);
    m = trans' * alpha(:,t-1);
    alpha(:,t) = m(:) .* obslik(:,t);
    [alpha(:,t), scale(t)] = normalise(alpha(:,t));
    %assert(approxeq(sum(alpha(:,t)),1))
end
loglik = sum(log(scale));

beta = zeros(Q,T);
t=T;
beta(:,T) = ones(Q,1);
gamma(:,T) = alpha(:,T);

for t=T-1:-1:1
    trans = reshape(transmat(:,:,t),[Q,Q]);
    b = beta(:,t+1) .* obslik(:,t+1);
    beta(:,t) = normalise(trans * b);
    gamma(:,t) = normalise(alpha(:,t) .* beta(:,t));
end