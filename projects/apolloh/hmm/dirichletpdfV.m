function p = dirichletpdfV(x, alpha)
%DIRICHLETPDF Dirichlet probability density function.
%   p = dirichletpdf(x, alpha) returns the probability of vector
%   x under the Dirichlet distribution with parameter vector
%   alpha.
%   
%
%   Author: David Ross
%   Modified by Sohrab Shah for vectorized version

%-------------------------------------------------
% Check the input
%-------------------------------------------------
error(nargchk(2,2,nargin));

% % enusre alpha is a vector
% if min(size(alpha)) ~= 1 | ndims(alpha) > 2 | length(alpha) == 1
%     error('alpha must be a vector');
% end

% ensure x is is a vector of the same size as alpha
if any(size(x) ~= size(alpha))
    error('x and alpha must be the same size');
end


%-------------------------------------------------
% Main
%-------------------------------------------------
% if any(x < 0,)
%     p = 0;
% elseif ~approxeq(sum(x),1,1e-3)
%     disp(['dirichletpdf warning: sum(x)~=1, but this may be ' ...
%         'due to numerical issues']);
%     disp(x)
%     disp(sum(x))
%     p = 0;
% else
    z = gammaln(sum(alpha)) - sum(gammaln(alpha)); 
    %z = exp(z);
    p = exp(z + sum((alpha-1).*log(x)));
    
    %.^(alpha-1));
    
%end