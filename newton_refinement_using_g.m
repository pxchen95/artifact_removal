function [w, d, num_iter, B, A, alp] = ...
    newton_refinement_using_g(w0, d0, maxiter, S, t, fs, K, tol)
% Solves the following least squares problem using Newton's descent: 
%           min_{w,d} min_{amp_0,amp_{i,c},amp_{i,s}} |S - A|^2,
% where 
% A(t) = amp_0 + sum_i^K amp_{i,c}*cos(2*pi*w*i*t) + amp_{i,s}*sin(2*pi*w*i*t)).
% In other words, we refine the frequency estimate w and phase shift
% estimates d while fitting A using harmonic regression. We remove the
% artifact as B = S - A.
%
% INPUTS: 
%     w0:         scalar, initial guess for frequency
%     d0:         1 x n vector, d(i) = initial guess for phase shift i
%     maxiter:    scalar, max # of iterations for Newton
%     S:          1 x (n+1) cell array, S{i} = 1 x N_i, samples of segment i
%     t:          1 x (n+1) cell array, t{i} = 1 x N_i, times that segment i is sampled at
%     fs:         scalar, sampling rate
%     K:          scalar, # of harmonics to fit
%     tol:        scalar, tolerance for stopping criteria
%
% OUTPUTS: 
%     w:        scalar, frequency
%     d:        1 x n vector, d(i) = phase shift i
%     num_iter: scalar, # of iterates used for Newton's ascent
%     B:        1 x (n+1) cell array, B{i} = 1 x N_i, recovered signal samples in segment i
%     A:        1 x (n+1) cell array, A{i} = 1 x N_i, reconstructed artifact samples in segment i
%     alp:      2*K+1 x 1 vector, 
%          alp(1)   = amp_0
%          alp(i)   = amplitude of cos(2*pi*K*i*t), i = 2, ..., K+1
%          alp(K+i) = amplitude of sin(2*pi*K*i*t), i = 2, ..., K+1

numSegments = length(S);
p0 = zeros(numSegments,1); % initialization
p0(1) = w0;
p0(2:end) = d0';
[~, ~, ~, grad_g, hess_g, g] = remove_artifact_ver_g(S, t, fs, K, p0(1), p0(2:end)');

% for debugging
p_save = zeros(numSegments, maxiter);
g_save = zeros(1,maxiter);
grad_g_save = zeros(numSegments, maxiter);
alp_1_save = zeros(1,maxiter);

p_save(:,1) = p0; 
g_save(1) = g;
grad_g_save(:,1) = grad_g;

for i = 1:maxiter
    dk = hess_g\grad_g;  % -descent direction
%     [alp_1, ~] = backtracking_linesearch_for_g(S, t, fs, K, p0(1), p0(2:end)', 100, dk);
    alp_1 = 1;           % step size
    p1 = p0 - alp_1*dk;  % Newton's descent
    
    err = norm(p1 - p0); % difference between iterates
    
    p0 = p1; % update iterate
    [B, A, alp, grad_g, hess_g, g] = remove_artifact_ver_g(S, t, fs, K, p0(1), p0(2:end)');

    % for debugging
    p_save(:,i+1) = p1;
    grad_g_save(:,i+1) = grad_g;
    g_save(i+1) = g;
    alp_1_save(i) = alp_1;
   
    % stopping criteria: found critical pt or insufficient change b/t iterates
    if norm(grad_g) < tol || err < eps
        break
    end
end

% outputs
w = p1(1);
d = p1(2:end)';
num_iter = i;

if num_iter == N_maxiter
    disp('WARNING: Newton refinement did not converge')
    fprintf('Stopping criteria values: gradient norm = %e, iterate diff = %e\n', ...
        norm(grad_g), err)
end

end