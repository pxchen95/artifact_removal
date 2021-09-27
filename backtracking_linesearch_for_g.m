function [alp_1, num_iter] = backtracking_linesearch_for_g(S, t, fs, K, p0, d0, maxiter, dk)
% Finds an appropriate step size by finding alp_1 that sufficiently 
% decreases g(x0 - alp_1*dk) relative to E(x0) (where x0 = [p0, d0])
% using backtracking linesearch
%
% INPUTS:
%     S:         1 x (n+1) cell array, S{i} = 1 x N_i, observed signal samples in segment i
%     t:         1 x (n+1) cell array, t{i} = 1 x N_i, unshifted times that segment i is sampled at, \in [0, T_i]
%     fs:        scalar, sample rate
%     K:         scalar, number of harmonics
%     p0:        scalar, current frequency iterate 
%     d0:        1 x n vector, current phase shift iterates
%     maxiter:   scalar, max # of iterations     
%     dk:        (n+1) x 1 vector, descent direction
%
% OUTPUTS:
%     alp_1:     scalar, step size
%     num_iter:  scalar, number of iterations

% linesearch control parameters
tau = 0.5;         % shrink factor for step size
c = 0.8;
[~, ~, ~, grad_g] = remove_artifact(S, t, fs, K, p0, d0);
m = -grad_g' * dk; % "local slope" of g(x0 - [.]*dk) along dk
param_t = -c*m;    % parameter for Armijo–Goldstein condition
[~, ~, ~, ~, ~, g_x] = remove_artifact(S, t, fs, K, p0, d0);

alp_1 = min(norm(dk),1); % initialize stepsize

% backtracking linesearch
for i = 1:maxiter
    [~, ~, ~, ~, ~, g_new] = remove_artifact(S, t, fs, K, p0 + alp_1*dk(1), d0 + alp_1*dk(2:end)');
    
    % Armijo–Goldstein condition || step size is too small
    if g_x - g_new > alp_1*param_t || alp_1 < 1e-8
        break
    end
    
    alp_1 = tau*alp_1; % shrink step size by tau 
end

num_iter = i;

if num_iter == maxiter
    disp('WARNING: backtracking linesearch for g did not converge')
    fprintf('Stopping criteria values: g_x - g_new = %e >? alp_1*param_t = %e; alp_1 = %f\n', g_x - g_new, alp_1*param_t, alp1)
end

end