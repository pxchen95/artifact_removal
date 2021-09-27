function [w, d, num_iter] = newton_ascent(w0, d0, N_maxiter, LS_maxiter, S, t, fs, tol)
% Find frequency w and phase shifts d that maximize E(w,d) using Newton's 
% ascent
% 
% INPUTS: 
%     w0:         scalar, initial guess for frequency
%     d0:         1 x n vector, d(i) = initial guess for phase shift i
%     N_maxiter:  scalar, max # of iterates for Newtons's ascent
%     LS_maxiter: scalar, max # of iterates for linesearch
%     S:          1 x (n+1) cell array, S{i} = 1 x N_i, samples of segment i
%     t:          1 x (n+1) cell array, t{i} = 1 x N_i, times that segment i is sampled at
%     fs:         scalar, sampling rate
%     tol:        scalar, tolerance for stopping criteria
%
% OUTPUTS: 
%     w:        scalar, frequency
%     d:        1 x n vector, d(i) = phase shift i
%     num_iter: scalar, # of iterates used for Newton's ascent

if ~iscell(S)
    disp('ERROR: Input S must be a cell array')
    return
end

if ~iscell(t)
    disp('ERROR: Input t must be a cell array')
    return
end

x0 = [w0; d0']; % initialize
[~, grad_E, hess_E] = calc_E_multtshifts(x0(1), x0(2:end)', S, t, fs);

for i = 1:N_maxiter
    % modified cholesky to ensure we use a negative definite matrix
    R = mod_cholesky(hess_E);
    dk = R\grad_E; % direction of ascent
    
    % backtracking linesearch to find step size
    [eta1, ~] = backtracking_linesearch(x0, dk, LS_maxiter, S, t, fs);
    
    % Newton's ascent
    x1 = x0 - eta1*dk; 
    
    % difference between iterates
    err = sum((x1 - x0).^2);
    
    x0 = x1; % update iterate
    [~, grad_E, hess_E] = calc_E_multtshifts(x0(1), x0(2:end)', S, t, fs);
    
    % stopping criteria: found critical pt and insufficient change b/t iterates
    if (sum(grad_E.^2) < tol && err < tol) || err < 1e-16
%         % for debugging
%         disp('grad norm')
%         sum(grad_E.^2)
%         disp('err')
%         err
        break
    end
end

% outputs
num_iter = i;
w = x1(1);
d = x1(2:end)';

if num_iter == N_maxiter
    disp('WARNING: Newton ascent did not converge')
    fprintf('Stopping criteria values: gradient norm squared = %e, iterate diff = %e\n', ...
        sum(grad_E.^2), err)
end

end