function [f, df] = grad_uv(uv, C2, yDims, R, lbda, include_mu)
%
% [f, df] = grad_uv(uv, C2, yDims, R, lbda)  
%
% Description: Gradient computation for low-rank alpha factorization 
%              optimization. This function is called by minimize.m.
%
% Arguments:
%
%     uv         -- variable with respect to which optimization is performed,
%                   where uv = [U(:); V(:); mu_u, mu_v]
%     C2         -- (numGroups x xDim) array; second moments of C
%     yDims      -- (1 x numGroups) array; dimensionalities of each
%                   observed group
%     R          -- int; rank of factorization
%     lbda       -- float; regularization parameter
%     include_mu -- logical; set true to include mean parameters in
%                   low-rank alpha model
%
% Outputs:
%
%     f          -- value of lower bound at uv
%     df         -- gradient at uv; same dimensions as uv    
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     01 Jul 2022 -- Initial full revision.  
%     09 Jul 2022 -- Added R = 0 compatibility.
%     12 Jul 2022 -- Added option to remove mean parameters from low-rank
%                    alpha model.

[numGroups, xDim] = size(C2);
lenU = numGroups*R;
lenV = xDim*R;

% Unpack uv
logalpha = zeros(numGroups,xDim);
if R > 0
    U = reshape(uv(1:lenU), [numGroups R]);
    V = reshape(uv(lenU+1:lenU+lenV), [xDim R]);
    logalpha = logalpha + U * V';
end

if include_mu
    mu_u = uv(lenU+lenV+1:lenU+lenV+numGroups);
    mu_v = uv(lenU+lenV+numGroups+1:end);
    logalpha = logalpha + mu_u*ones(1,xDim) + ones(numGroups,1)*mu_v';
end

f = trace(logalpha * ones(xDim,1) * yDims) ...
    - trace(C2 .* exp(logalpha) * ones(xDim,numGroups));
if R > 0
    f = f - lbda * trace(U'*U + V'*V);
end
    
A = yDims'*ones(1,xDim) - C2.*exp(logalpha); % Intermediate value

df = [];
if R > 0
    dU = A*V - 2*lbda*U;
    dV = A'*U - 2*lbda*V;
    df = [df; dU(:); dV(:)];
end
if include_mu
    dmu_u = A*ones(xDim,1);
    dmu_v = A'*ones(numGroups,1);
    df = [df; dmu_u; dmu_v];
end

% Account for final constants and the fact that we're minimizing
f = -(f/2);
df = -(df/2);
