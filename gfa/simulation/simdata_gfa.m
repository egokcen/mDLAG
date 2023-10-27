function out = simdata_gfa(N, yDims, xDim, hyperparams, snr, varargin)
% 
% out = simdata_gfa(N, yDims, xDim, hyperparams, snr, ...)
%
% Description: Generate simulated data according to a group factor analysis
%              model.
%
% Arguments:
%
%     Required:
%
%     N        -- int; number of data points to generate
%     yDims    -- (1 x M) array; List of dimensionalities of
%                 observed data, [y1Dim, y2Dim, ...]
%     xDim     -- int; Dimensionality of latents, X
%     hyperparams -- structure with the following fields:
%                      beta  -- positive float; precision of mean parameter
%                               generative model (Gaussian)
%                      a_phi -- positive float; 'a' shape parameter of
%                               observation precision (phi) generative 
%                               model (Gamma with mean a/b)
%                      b_phi -- positive float; 'b' scale parameter of
%                               observation precision (phi) generative
%                               model (Gamma with mean a/b)
%                    if R == 'full':
%                      a_alpha -- (M x xDim) array; 'a' shape parameter of 
%                                 alpha parameter generative model 
%                                 (Gamma with mean a/b)
%                      b_alpha -- (M x xDim) array; 'b' scale parameter of
%                                 alpha parameter generative model 
%                                 (Gamma with mean a/b)
%                    if R != 'full':
%                      uv_lbda -- positive float; precision for low-rank
%                                 alpha parameters
%     snr      -- (1 x numGroups) array; List of signal-to-noise ratios,
%                 defined as trace(CC') / sum(1./phi)
%     Optional:
%
%     R          -- string or int; 'full' for full-rank alpha model, a
%                   positive integer giving the rank of alpha otherwise.
%                   (default: 'full')
%     include_mu -- logical; set true to include mean parameters in
%                   low-rank alpha model. (default: false)
%
% Outputs:
%     out -- output structure with the following fields:
%         Ys     -- (1 x M) cell array; list of data matrices 
%                   {(y1Dim x N), (y2Dim x N), ...}
%         X      -- (xDim x N) array; latent data
%         phis   -- (1 x M) cell array; List of observation 
%                   precision parameters {(y1Dim x 1), (y2Dim x 1), ...}
%         ds     -- (1 x M) cell array; List of data means
%                   {(y1Dim x 1), (y2Dim x 1), ...}
%         alphas -- (M x xDim) array; alpha parameter values
%         Cs     -- (1 x M) cell array; List of factor loadings 
%                   {(y1Dim x xDim), (y2Dim x xDim), ...}
%       if R != 'full':
%         if R > 0
%           U     -- (M x R) array; low-rank factorization of alpha
%           V     -- (xDim x R) array; low-rank factorization of alpha
%         if include_mu
%           mu_u  -- (M x 1) array; mean offset in low-rank factorization
%           mu_v  -- (xDim x 1) array; mean offset in low-rank factorization
%
% Author: Evren Gokcen
%
% Revision History:
%     06 Jun 2022 -- Initial full revision.
%     03 Jul 2022 -- Added control over signal-to-noise ratio.
%     09 Jul 2022 -- Added low-rank option.
%     12 Jul 2022 -- Added option to remove mean parameters from low-rank
%                    alpha model.

R = 'full';
include_mu = false;
assignopts(who, varargin);

M = length(yDims);

if ~isequal(R, 'full')
    
    if R <= 0 && ~include_mu
        fprintf('Error: If R = 0, then mu must be included in the low-rank alpha model.\n');
        out = [];
        return
    end
    
    if R >= min([xDim M])
        fprintf('Specified R corresponds to full rank model.\n');
        R = 'full';
    end
end

% Initialize output structure
out = struct('Ys', [], ...
             'X', nan(xDim,N), ...
             'phis', [], ...
             'ds', [], ...
             'alphas', nan(M,xDim), ...
             'Cs', []);
out.Ys = cell(1,M);
out.phis = cell(1,M);
out.ds = cell(1,M);
out.Cs = cell(1,M);
if ~isequal(R, 'full')
    if R > 0
        out.U = nan(M,R);
        out.V = nan(xDim,R);
    end
    if include_mu
        out.mu_u = nan(M,1);
        out.mu_v = nan(xDim,1);
    end
end


% Generate latent data
out.X = randn(xDim, N);

% Alpha model, full or low rank
if isequal(R, 'full')
    for m = 1:M
        for xIdx = 1:xDim
            out.alphas(m,xIdx) = gamrnd(hyperparams.a_alpha(m,xIdx), 1./hyperparams.b_alpha(m,xIdx));
        end
    end
else
    % Generate low rank alpha parameters
    out.alphas = zeros(M,xDim);
    if R > 0
        out.U = hyperparams.uv_lbda^(-1).*randn(M,R);
        out.V = hyperparams.uv_lbda^(-1).*randn(xDim,R);
        out.alphas = out.alphas + out.U * out.V';
    end
    if include_mu
        out.mu_u = hyperparams.uv_lbda^(-1).*randn(M,1);
        out.mu_v = hyperparams.uv_lbda^(-1).*randn(xDim,1);
        out.alphas = out.alphas + out.mu_u*ones(1,xDim) + ones(M,1)*out.mu_v';
    end
    out.alphas = exp(out.alphas);
end

% Remaining parameters
for m = 1:M
    
    % Generate parameters
    out.ds{m} = mvnrnd(zeros(1,yDims(m)), hyperparams.beta^(-1).*eye(yDims(m)), 1)';
    out.phis{m} = gamrnd(hyperparams.a_phi, 1./hyperparams.b_phi, yDims(m), 1);
    out.Cs{m} = nan(yDims(m),xDim);
    for xIdx = 1:xDim
        out.Cs{m}(:,xIdx) = mvnrnd(zeros(1,yDims(m)), out.alphas(m,xIdx)^(-1).*eye(yDims(m)), 1)';
    end
    
    % Enforce the desired signal-to-noise ratios
    varCC = trace(out.Cs{m} * out.Cs{m}');
    varNoise_desired = varCC / snr(m);
    varNoise_current = sum(out.phis{m}.^(-1));
    out.phis{m} = out.phis{m} .* (varNoise_current / varNoise_desired);
    
    % Generate observations
    ns = mvnrnd(zeros(1,yDims(m)), diag(out.phis{m}.^(-1)), N)';
    out.Ys{m} = out.Cs{m} * out.X + repmat(out.ds{m},1,N) + ns;
end