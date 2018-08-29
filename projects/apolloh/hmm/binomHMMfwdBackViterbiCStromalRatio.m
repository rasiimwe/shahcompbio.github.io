function [G,sReturn,rtReturn,muReturn,piReturn,ll,rho,segs] = binomHMMfwdBackViterbiCStromalRatio(ref,depth,chr,posn,cn,r_0,rn,estimateS,s_0,alphaHyperR,betaHyperR,alphaHyperS,betaHyperS,pi_0,kappa,txnExpLen,ZS,CNS,maxiter)
%% FUNCTION [G,sReturn,rtReturn,muReturn,piReturn,ll,rho,segs] = binomHMMfwdBackViterbiCStromalRatio(ref,depth,chr,posn,cn,r_0,rn,s_0,alphaHyperS,betaHyperS,pi_0,kappa,txnExpLen,ZS,CNS,maxiter)
%
% author: Gavin Ha <gha@bccrc.ca>
%         Dept of Molecular Oncolgy/Centre for Translational and Applied Genomics
%         British Columbia Cancer Agency
%         University of British Columbia
% date  : March 7, 2011

K = length(r_0); % prior probability of success
N = length(ref); % number of data points
py = zeros(K,N); % local evidence
mus = zeros(K,maxiter); % state Binomial means
pi = zeros(K,maxiter); % initial state distribution 
s = zeros(1,maxiter); % stromal parameter
rt = zeros(K,maxiter); % tumour reference allelic ratio (1 for each state)
converged = 0; % flag for convergence
G = zeros(N,1); %latent variable
loglik = zeros(1,maxiter); %log likelihood

%% SET UP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up the chromosome indicies
% make cell array of chromosome indicies
chrs = unique(chr);
chrsI = cell(1,length(chrs));

piZ = cell(1,length(chrs)); % need 1 initial dist per chromosome
% initialise the chromosome index and the init state distributions 
for c = 1:length(chrs)
    chrsI{c} = find(chr == chrs(c));
    piZ{c} = ones(1,K)/K;
end

%% INITIALISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 1;
s(i) = s_0;
rt(:,i) = r_0;
mus(:,i) = stromalTwoComponentMixture(r_0,rn,s_0);
pi(:,1) = pi_0;

% CY matrix; KxN
%disp('binomHMMfwdBackViterbiC: Initializing genotype switch matrix');
%CY = zeros(K,N);
%CY = genotypeSwitchMatrixC(CY,cn,CNS);

%find positions with depth higher than 1000
highDepthInd = depth > 1000;

% calculate the likelihood conditional on k
for k=1:K
    py(k,:)=binomialpdf(ref,depth,mus(k,i));% + 2.225e-307;    
end
py(:,highDepthInd) = py(:,highDepthInd) + 2.225e-307; %add pseudocounts to high depth positions to avoid underflow
%py = py .* CY;
%py = normalise(py .* CY,1);

% initialise transition matrix to the prior:
disp('binomHMMfwdBackViterbiC: Initializing Non-stationary transition matrix');
A = zeros(K,K);

%% Expectation Maximization via Forwards-Backwards
loglik(i) = -Inf;
rho = zeros(K,N);    % marginal p(Z_t|Y)
%xi  = zeros(K,K,N);  % marginal p(Z_t,Z_{t-1}|Y); used for estimating A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (~converged && (i < maxiter))
    disp(['binomHMMfwdBackViterbiC: EM iteration:', int2str(i), ' loglik: ',num2str(loglik(i))]);
    i = i+1;
    % E-step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for c = 1:length(chrsI)
        piZ{c} = pi(:,i-1)';
        I = chrsI{c};
        % function [gamma,alpha,beta,loglik] = nonStationaryFwdBack(init_state_distrib,transmat,obslik)
        [rho(:,I),a,b,ll] = fwd_backC_cn(piZ{c}, A, py(:,I),cn(I),CNS,ZS,posn(I),txnExpLen);       
        loglik(i) = loglik(i)+ll;
    end
    %Zcounts = reshape(sum(xi,3),[K,K]);
    
%    disp(Zcounts);
    % M-step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update the noise hyperparams
    [s(i), rt(:,i), pi(:,i)] = ...
        estimateBinomNoiseStromalRatioParamsMap(ref',depth',rho,estimateS,rt(:,i-1),rn,alphaHyperR,betaHyperR,alphaHyperS,betaHyperS,kappa);

    % re-calculate the likelihood conditional on k
    mus(:,i) = stromalTwoComponentMixture(rt(:,i),rn,s(i));
    for k=1:K
        py(k,:)=binomialpdf(ref,depth,mus(k,i));% + 2.225e-307;
    end
    py(:,highDepthInd) = py(:,highDepthInd) + 2.225e-307; %add pseudocounts to high depth positions to avoid underflow
    %py = py .* CY;
    %py = normalise(py .* CY,1);
 
    % compute log-likelihood and check convergence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    for k=1:K
        priorR(k) = betapdflog(rt(k,i),alphaHyperR(k),betaHyperR(k));        
    end
    if estimateS==1 %if specified to estimate stromal contamination
        priorS = betapdflog(s(i),alphaHyperS,betaHyperS);    
    else
        priorS = 0;
    end
    priorPi = dirichletpdflog(pi(:,i)',kappa);
    
    disp(['priorS ',num2str(priorS)])
    disp(['priorR ',num2str(sum(priorR))])
    disp(['priorPi ',num2str(priorPi)])
    disp(['loglik ',num2str(loglik(i))])
    loglik(i) = loglik(i) + priorS + sum(priorR) + priorPi;
    if approxeq(loglik(i),loglik(i-1),1)% || (loglik(i) < loglik(i-1))
        converged = 1;
    elseif (loglik(i) < loglik(i-1))
        converged = 1;
        i = i-1;
        % recompute Binomial parameter using converged iteration i-1 before
        % exiting EM       
        for k=1:K
            py(k,:)=binomialpdf(ref,depth,mus(k,i));% + eps;%2.225e-307;
        end
        py(:,highDepthInd) = py(:,highDepthInd) + 2.225e-307; %add pseudocounts to high depth positions to avoid underflow
    end    
end
%% Perform one last round of E-step to get latest responsibilities
% E-step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for c = 1:length(chrsI)
    piZ{c} = pi(:,i-1)';
    I = chrsI{c};
    % function [gamma,alpha,beta,loglik] = nonStationaryFwdBack(init_state_distrib,transmat,obslik)
    [rho(:,I),a,b,ll] = fwd_backC_cn(piZ{c}, A, py(:,I),cn(I),CNS,ZS,posn(I),txnExpLen);       
    loglik(i) = loglik(i)+ll;
end

%% Viterbi
disp(['Using parameters from EM iteration ',num2str(i)]);        
segs = cell(1,length(chrs));
% do the final viterbi pass
for c = 1:length(chrsI)
    piZ{c} = pi(:,i)';
    I = chrsI{c};
    %function [path loglik] = nonStationaryViterbiPath(prior, transmat,obslik)
    [G(I) ll segs{c}] = viterbi_pathC_cn(log(piZ{c}), A, log(py(:,I)), cn(I), CNS, ZS, posn(I), txnExpLen);
end

sReturn = s(1:i);
rtReturn = rt(:,1:i);
piReturn = pi(:,1:i);
muReturn  = mus(:,1:i);
ll = loglik(1:i);
