function A = initializeNSTxnKKT(ZS,CNS,chr,posn,cn,txnExpLen)
%% FUNCTION A = initializeNSTxnKKT(ZS,CNS,chr,posn,cn,txnExpLen)
%
% author: Gavin Ha <gha@bccrc.ca>
%         Dept of Molecular Oncolgy/Centre for Translational and Applied Genomics
%         British Columbia Cancer Agency
%         University of British Columbia
% date  : April 12, 2011

K = length(ZS);
N = length(chr);
A = zeros(K,K,N);
%C = repmat(cnUncertainty,[K,K]);
for t = 2:N  
    if (chr(t)==chr(t-1))
        d = posn(t) - posn(t-1) + 1;
    else
        d = posn(t);
    end
    txnRho = 1 - distanceTransitionFunction(d,txnExpLen);
    prevC = CNS(cn(t-1)):(CNS(cn(t-1)+1)-1);
    curC = CNS(cn(t)):(CNS(cn(t)+1)-1);
    C = zeros(K,K);

    % now set 1's to allowable transitions
    for ci = prevC
        for cj = curC
           C(ci,cj) = 1; 
        end
    end
       
    
    for ai = 1:K
        %rowSum = 0;
        for aj = 1:K
            if (ZS(ai)==ZS(aj))
                A(ai,aj,t) = txnRho;                
            else
                A(ai,aj,t) = (1-txnRho)/(length(curC));
            end
            %rowSum = rowSum + A(ai,aj,t);
        end
%         for aj = curC
%             A(ai,aj,t) = A(ai,aj,t) / rowSum;
%         end
    end   
    A(:,:,t) = A(:,:,t) .* C;
    A(:,:,t) = normalise(A(:,:,t),2);    
end