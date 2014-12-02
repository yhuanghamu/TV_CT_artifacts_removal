function [X,info] = TVdenoise(B,delta,eps_rel)
%TVDENOISE  Total variation image denoising
% 
% X = TVdenoise(B,delta)
% [X,info] = TVdenoise(B,delta,eps_rel)
%
% This function solves the TV denoising problem
%
%    min  TV(X)  subject to   || X - B ||_F <= delta
%
% where B is a noisy image, X is the reconstruction, and delta is an
% upper bound for the residual norm.  The TV function is the 1-norm
% of the gradient magnitude, computed via neighbor pixel differences.
% At the image borders, we imposed reflexive boundary conditions for
% the gradient computations.
%
% The parameter delta should be of the same size as the norm of the
% image noise.  If the image is m-times-n, and sigma is the standard
% deviation of the image noise in a pixel, then we recommend to use
% delta = tau*sqrt(m*n)*sigma, where tau is slightly smaller than one,
% say, tau = 0.85.
%
% The function returns an epsilon-optimal solution X, meaning that
% if X* is the exact solution, then our solution X satisfies
%
%     TV(X) - TV(X*) <= epsilon = max(B(:))*m*n*eps_rel,
%
% where eps_rel is a specified relative accuracy (default eps_rel = 1e-3).
% 
% The solution status is returned in the stuct info, with info.STATUS
% having one of the settings
%  'EPSILON-OPTIMAL-SOLUTION': X is an epsilon-optimal solution
%  'NOT-EPSILON-OPTIMAL-SOLUTION': X is not an epsilon-optimal solution
%  'MAXIMUM-NUMBER-OF-ITERATIONS-EXCEEDED': X is not an epsilon-optimal
%     solution when the maximum number of iterations was reached. 
% Other fields of info:
%  info.NDENOISE      Upper bound for the number of iterations.
%  info.ITERATIONS_K  The number of iterations used.	
%  info.EPS_REL_K     The relative accuracy reached.
%  info.TIME          Time in seconds for the program to run. 
%
% See also: TVdeblur, TVinpaint.

% J. Dahl^1, P.C. Hansen^2, S.H. Jensen^1 & T.L. Jensen^1
% CSI project: (1) Aalborg University, (2)Technical University of Denmark
% April 28, 2009.

% Check input parameters.
if nargin < 2	
	error('Too few input parameters');
elseif nargin == 2
	eps_rel = 1e-3;
end
tic

% The special case where the solution is the constant image.
alpha = sum(B(:))/numel(B);
X = alpha*ones(size(B));
mdelta = norm(X-B,'fro');
if  mdelta < delta
    % The constant image is the solution.
    info = info_type_denoise(1,0,0,0,toc);
    return;
elseif mdelta < 1.1*delta
    % The constant image is almost a solution.
    warning('Convergence problems may arise')
end

% Set parameters for denoising algorithm.
R = max(B(:));
mn = numel(B);
epsilon = R*mn*eps_rel;% duality gap
mu = epsilon/mn;% regularization arameter
Lmu = 8/mu ;% Lipschitz continuous derivatives
N = int32( ceil(2*sqrt(8*mn)*delta/epsilon) );% iteration upper bound

% Compute TV solution via C function.
[X,k,epsilon_k] = tv_denoise(B,delta,epsilon,Lmu,mu,N,0);%k:current iteration number; epsilon_k:current epsilon.

% Set info, if required.
if nargout == 2
	if k >= N
	    info = info_type_denoise(3,N,k,epsilon_k/(R*mn),toc);
    elseif epsilon_k > epsilon
		info = info_type_denoise(2,N,k,epsilon_k/(R*mn),toc);
	else
		info = info_type_denoise(1,N,k,epsilon_k/(R*mn),toc);	
    end
end