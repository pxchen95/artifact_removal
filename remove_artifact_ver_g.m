function [B, A, alp, grad_g, hess_g, g] = remove_artifact_ver_g(S, t, fs, K, freq, delta_hat)
% (1) Reconstructs periodic artifact A with fundamental frequency freq 
% using K harmonics, i.e. A(t) = sum_i^K amp_{i,c}*cos(2*pi*w*i*t) +
% amp_{i,s}*sin(2*pi*w*i*t)), and removes it from S, i.e. B = S - A.
% (2) Computes objective function g in least squares fitting:
%        g(w,d) = min_{amp_0,amp{i,c},amp_{i,s}} |S - A|^2,
% and its gradient grad_g and hessian hess_g at (w,d) = (freq, delta_hat).
% 
% INPUTS: 
%     S:         1 x (n+1) cell array, S{i} = 1 x N_i, observed signal samples in segment i
%     t:         1 x (n+1) cell array, t{i} = 1 x N_i, unshifted times that segment i is sampled at, \in [0, T_i]
%     fs:        scalar, sample rate
%     K:         scalar, number of harmonics
%     freq:      scalar, true/estimated frequency
%     delta_hat: 1 x n vector, true/estimated phase shifts
%
% OUTPUTS:
%     B:      1 x (n+1) cell array, B{i} = 1 x N_i, recovered signal samples in segment i
%     A:      1 x (n+1) cell array, A{i} = 1 x N_i, reconstructed artifact samples in segment i
%     alp:    2*K+1 x 1 vector, 
%          alp(1)   = amp_0
%          alp(i)   = amplitude of cos(2*pi*K*i*t), i = 2, ..., K+1
%          alp(K+i) = amplitude of sin(2*pi*K*i*t), i = 2, ..., K+1
%     grad_g: (n+1) x 1 vector, gradient of g at (freq, delta_hat)
%     hess_g: (n+1) x (n+1) matrix, hessian of g at (freq, delta_hat)
%     g:      scalar, value of least squares function g at (freq, delta_hat)

% CONVERT CELL ARRAYS TO VECTORS
numSegments = length(S); % number of segments
S_vec = []; t_vec = []; t_vec_no_d = [];
N = zeros(1, numSegments);
sumT = 0;
for i = 1:numSegments
    S_vec = [S_vec, S{i}];  % observed signal
    
    if i == 1
        sum_shift = 0;
    else
        sum_shift = delta_hat(i-1)/freq; % sum of time shifts
    end
    
    N_i = length(S{i});        % number of samples in segment i
    T_i = N_i/fs;              % length of segment i in time
    tshift = sumT + sum_shift; % segment i starts sampling at time tshift
    
    N(i) = N_i; % length of segment i in samples
    
    t_vec = [t_vec, t{i}+tshift]; % sample times (shifted)
    t_vec_no_d = [t_vec_no_d, t{i}+sumT]; % used to calculate grad_g, hess_g
    
    sumT = sumT + T_i; % update endtime sum
end

% COMPUTE AMPLITUDES USING LEAST SQUARES IN TIME
cos_t = @(k) cos(2*pi*freq*k*t_vec); 
sin_t = @(k) sin(2*pi*freq*k*t_vec);

% construct matrix
A = zeros(2*K,2*K);
for i = 1:2*K
    for j = 1:i
        if i <= K && j <= K
            A(i,j) = cos_t(i)*cos_t(j)';
        elseif i > K && j > K
            A(i,j) = sin_t(i-K)*sin_t(j-K)';
        elseif i <= K && j > K
            A(i,j) = cos_t(i)*sin_t(j-K)';
        else
            A(i,j) = sin_t(i-K)*cos_t(j)';
        end
        
        A(j,i) = A(i,j);
    end
end
A_new = zeros(2*K+1,2*K+1);
A_new(2:end, 2:end) = A;
A_new(1,1) = length(t_vec);

A_row1 = zeros(1,2*K);
for i = 1:2*K
    if i <= K
        A_row1(i) = sum(cos_t(i));
    else
        A_row1(i) = sum(sin_t(i-K));
    end
    A_new(1,i+1) = A_row1(i);
    A_new(i+1,1) = A_new(1,i+1);
end
A = A_new;

% right-hand side
b = zeros(2*K, 1);
for i = 1:2*K
    if i <= K
        b(i) = S_vec*cos_t(i)';
    else
        b(i) = S_vec*sin_t(i-K)';
    end
end
b = [sum(S_vec); b];

% amplitudes of harmonics
alp = A\b;

% RECONSTRUCT ARTIFACT
alp_0 = alp(1);
alp = alp(2:end);
A_vec = zeros(size(S_vec)); 
g_coeff = zeros(size(t_vec));
g_coeff2 = g_coeff; 
for i = 1:K
    A_vec = A_vec + alp(i)*cos_t(i) + alp(K+i)*sin_t(i); 
    g_coeff = g_coeff + i*(alp(i)*sin_t(i) - alp(K+i)*cos_t(i));
    g_coeff2 = g_coeff2 + i^2*(alp(i)*cos_t(i) + alp(K+i)*sin_t(i));
end
g_coeff = g_coeff * (2*pi);
g_coeff2 = g_coeff2 * (2*pi)^2;
A_vec = A_vec + alp_0;
alp = [alp_0; alp];

% REMOVE ARTIFACT
B_vec = S_vec - A_vec; 

% COMPUTE g, grad_g, hess_g
g = norm(B_vec)^2; 

grad_g = zeros(numSegments,1);
grad_g(1) = 2*sum(B_vec.*g_coeff.*t_vec_no_d);

hess_g = zeros(numSegments,numSegments);
hess_g(1,1) = 2*sum((g_coeff.*t_vec_no_d).^2 + B_vec .* g_coeff2 .* t_vec_no_d.^2);

% CONVERT VECTORS TO CELL ARRAYS
B = {}; A = {}; g_coeff_cell = {}; g_coeff2_cell = {}; t_vec_no_d_cell = {};
ind1 = 1; 
for i = 1:numSegments
    ind2 = ind1 + N(i) - 1;
    B{i} = B_vec(ind1:ind2);
    A{i} = A_vec(ind1:ind2);
    g_coeff_cell{i} = g_coeff(ind1:ind2);
    g_coeff2_cell{i} = g_coeff2(ind1:ind2);
    t_vec_no_d_cell{i} = t_vec_no_d(ind1:ind2);
    ind1 = ind2 + 1;
end

% COMPUTE grad_g, hess_g
for i = 2:numSegments
    grad_g(i) = 2*sum(B{i} .* g_coeff_cell{i});
    hess_g(1,i) = 2*sum((t_vec_no_d_cell{i} .* g_coeff_cell{i}.^2) ...
        + B{i} .* g_coeff2_cell{i} .* t_vec_no_d_cell{i});
    hess_g(i,1) = hess_g(1,i);
    hess_g(i,i) = 2*sum(g_coeff_cell{i}.^2 + B{i}.*g_coeff2_cell{i});
end

end