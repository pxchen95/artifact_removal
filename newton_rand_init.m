function [w, d, num_iter, t, w0_test, d0_test, E_save, freq_save, d_save, NA_num_iter_save] = ...
    newton_rand_init(w0, w0_range, numInit, N_maxiter, LS_maxiter, S, fs, tol)
% Run newton_ascent.m using (uniform) random initialization
% 
% INPUTS: 
%     w0:         scalar, initial guess for frequency
%     w0_range:   scalar, search for freq initializations in [w0 - w0_range, w0 + w0_range]
%     numInit:    scalar, # of random initializations to test (sugg. 10)
% NOTE: all inputs below are inputs for newton_ascent.m
%     N_maxiter:  scalar, max # of iterates for Newtons's ascent
%     LS_maxiter: scalar, max # of iterates for linesearch
%     S:          1 x (n+1) cell array, S{i} = 1 x N_i, samples of segment i
%     fs:         scalar, sampling rate
%     tol:        scalar, tolerance for stopping criteria
%
% OUTPUTS: 
%     w:         scalar, frequency
%     d:         1 x n vector, d(i) = phase shift i
%     num_iter:  scalar, # of iterates used for Newton's ascent
%     t:         1 x (n+1) cell array, t{i} = 1 x N_i, "unshifted" sample times in [0, (N_i-1)/fs]
% NOTE: all outputs below are for debugging purposes
%     w0_test:   1 x numInit vector, freq_save(i) = frequency initializations used
%     d0_test:   numInit x numShifts matrix, shift initializations used
%     E_save:    1 x numInit vector, E(i) = energy corresponding to initialization i 
%     freq_save: 1 x numInit vector, frequency estimates for each initialization
%     d_save:    numInit x numShifts matrix, shift estimates for each initialization
%     NA_num_iter_save: 1 x numInit vector, # of iterates of Newton's ascent used for each initialization

if ~iscell(S)
    disp('ERROR: Input S must be a cell array')
    return
end

numSegments = length(S);     % # of segments
numShifts = numSegments - 1; % # of shifts = # of segments - 1

% compute t = 1 x (n+1) cell array, t{i} = 1 x N_i, "unshifted" sample times in [0, (N_i-1)/fs]
t = {};
for i = 1:numSegments
    t{i} = 0:length(S{i})-1; t{i} = t{i}/fs;
end

% random initialization
w0_test = w0 + rand(1,numInit)*2*w0_range - w0_range; % pick random inits for freq
w0_test(1) = w0;                                      % force 1 init to be the initial guess
d0_test = rand(numInit,numShifts);                    % pick random inits for shifts

freq_save = zeros(1,numInit); % save outputs of each initialization
E_save = freq_save; 
NA_num_iter_save = freq_save;
d_save = zeros(numInit,numShifts);

for i = 1:numInit % run Newton's ascent for each initialization
%     i
    [freq_tmp, d_tmp, NA_num_iter] = newton_ascent(w0_test(i), d0_test(i,:), N_maxiter, LS_maxiter, S, t, fs, tol);
    
    % save energy corresponding to output
    E_tmp = calc_E_multtshifts(freq_tmp, d_tmp, S, t, fs);    
    
    % save outputs
    freq_save(i) = freq_tmp; d_save(i,:) = d_tmp; E_save(i) = E_tmp; NA_num_iter_save(i) = NA_num_iter;
end

% pick "best" initialization = initialization that maximizes energy
[~,I] = max(E_save);

% corresponding outputs of best initialization used
w = freq_save(I);
d = d_save(I,:);
num_iter = NA_num_iter_save(I);

end