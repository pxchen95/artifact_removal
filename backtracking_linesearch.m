function [eta1, num_iter] = backtracking_linesearch(x0, d0, N_maxiter, f, t, fs)
% Finds an appropriate step size by finding eta1 that sufficiently 
% increases E(x0 - eta1*d0) relative to E(x0) using backtracking linesearch
%
% INPUTS:
%     x0:        (n+1) x 1 vector, current iterate
%     d0:        (n+1) x 1 vector, ascent direction
%     N_maxiter: scalar, max # of iterations     
%     f:         1 x (n+1) cell array, f{i} = 1 x N_i, samples of segment i
%     t:         1 x (n+1) cell array, t{i} = 1 x N_i, times that segment i is sampled at
%     fs:        scalar, sampling rate
%
% OUTPUTS:
%     eta1:     scalar, step size
%     num_iter: scalar, number of iterations
    
% max candidate step size 
if length(d0) > 1
    max_eta = 1/norm(d0(2:end)); 
else % no time gaps
    max_eta = 1;
end
eta1 = max_eta; % initialize

% linesearch control parameters
tau = 0.5;  % shrink factor for step size
c = 0.5;
[E, grad_E, ~] = calc_E_multtshifts(x0(1), x0(2:end)', f, t, fs);
E = -E;             % min(-E) = max(E)
grad_E = -grad_E;   % grad(-E) = -grad(E)
m =  -grad_E' * d0; % "local slope" of E(x0 - [.]*d0) along d0
alp = -c*m; % parameter for Armijo–Goldstein condition

% backtracking linesearch
for i = 1:N_maxiter
    [E0, ~, ~] = calc_E_multtshifts(x0(1)-eta1*d0(1), (x0(2:end)-eta1*d0(2:end))', f, t, fs);
    E0 = -E0; % min(-E) = max(E)

    % Armijo–Goldstein condition
    if E - E0 >  eta1*alp
        break
    end

    eta1 = eta1*tau; % shrink step size by tau
end

num_iter = i;

if num_iter == N_maxiter
    disp('WARNING: backtracking linesearch did not converge')
    fprintf('Stopping criteria values: E - E0 = %f >? eta1*alp = %f\n', E - E0, eta1*alp)
end

% constrain eta1 in [0, 1/norm(d0(2:end))]
if eta1 > max_eta || eta1 <= 0
    eta1 = max_eta;
end

end