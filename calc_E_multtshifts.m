function [E, grad_E, hess_E] = calc_E_multtshifts(w, d, f, t, fs)
% Computes energy E, its gradient grad_E, and its Hessian hess_E at(w,d) 
% using trapezoidal rule to approximate the integrals
%
% INPUTS:
%   w: scalar, frequency
%   d: 1 x n vector, d(i) = hat delta_i, phase shifts
%   f: 1 x (n+1) cell array, f{i} = 1 x N_i, samples of segment i
%   t: 1 x (n+1) cell array, t{i} = 1 x N_i, times that segment i is sampled at
%   fs: scalar, sampling rate
%
% OUTPUTS:
%   E:      scalar, energy E at (w,d)
%   grad_E: (n+1) x 1 vector, gradient of E at (w,d)
%   hess_E: (n+1) x (n+1) matrix, Hessian of E at (w,d)
    
d = [0, d];
n = length(t) - 1; % number of time shifts
T = zeros(1,n+1);  % length in time of segment i, T(i) = T_i

% COMPUTE COS/SIN
cos_wave_shift = {};
sin_wave_shift = {};
for i = 1:n+1
    cos_wave_shift{i} = cos(2*pi*(w*(t{i}+sum(T))+d(i)));
    sin_wave_shift{i} = sin(2*pi*(w*(t{i}+sum(T))+d(i)));
    T(i) = length(t{i})/fs;
end

% COMPUTE COS/SIN * F_i
cos_wave_shift_times_fi = cos_wave_shift;
sin_wave_shift_times_fi = sin_wave_shift;
for i = 1:n+1
    cos_wave_shift_times_fi{i} = cos_wave_shift_times_fi{i}.*f{i};
    sin_wave_shift_times_fi{i} = sin_wave_shift_times_fi{i}.*f{i};
end

T_cumsum = [0, cumsum(T)]; % T_cumsum(i) = sum_{j=1}^{i-1} T(j)

%  COMPUTE COS/SIN * F_i * t (or t^2)
cos_wave_shift_times_fi_times_t = cos_wave_shift_times_fi;
sin_wave_shift_times_fi_times_t = sin_wave_shift_times_fi;
cos_wave_shift_times_fi_times_t2 = cos_wave_shift_times_fi;
sin_wave_shift_times_fi_times_t2 = sin_wave_shift_times_fi;
for i = 1:n+1
    cos_wave_shift_times_fi_times_t{i} = cos_wave_shift_times_fi_times_t{i}.*(t{i} + T_cumsum(i));
    sin_wave_shift_times_fi_times_t{i} = sin_wave_shift_times_fi_times_t{i}.*(t{i} + T_cumsum(i));

    cos_wave_shift_times_fi_times_t2{i} = cos_wave_shift_times_fi_times_t{i}.*(t{i} + T_cumsum(i));
    sin_wave_shift_times_fi_times_t2{i} = sin_wave_shift_times_fi_times_t{i}.*(t{i} + T_cumsum(i));
end

% COMPUTE INTEGRALS USING TRAPEZOIDAL RULE
int_cos_wave_shift_times_fi = zeros(1,n+1);
int_sin_wave_shift_times_fi = zeros(1,n+1);
int_cos_wave_shift_times_fi_times_t = zeros(1,n+1);
int_sin_wave_shift_times_fi_times_t = zeros(1,n+1);
int_cos_wave_shift_times_fi_times_t2 = zeros(1,n+1);
int_sin_wave_shift_times_fi_times_t2 = zeros(1,n+1);
for i = 1:n+1
    cos_tmp = cos_wave_shift_times_fi{i};
    cos_tmp(1) = cos_tmp(1)/2; cos_tmp(end) = cos_tmp(end)/2;
    int_cos_wave_shift_times_fi(i) = sum(cos_tmp)/fs;

    sin_tmp = sin_wave_shift_times_fi{i};
    sin_tmp(1) = sin_tmp(1)/2; sin_tmp(end) = sin_tmp(end)/2;
    int_sin_wave_shift_times_fi(i) = sum(sin_tmp)/fs;

    % times t
    cos_tmp = cos_wave_shift_times_fi_times_t{i};
    cos_tmp(1) = cos_tmp(1)/2; cos_tmp(end) = cos_tmp(end)/2;
    int_cos_wave_shift_times_fi_times_t(i) = sum(cos_tmp)/fs;

    sin_tmp = sin_wave_shift_times_fi_times_t{i};
    sin_tmp(1) = sin_tmp(1)/2; sin_tmp(end) = sin_tmp(end)/2;
    int_sin_wave_shift_times_fi_times_t(i) = sum(sin_tmp)/fs;

    % times t^2
    cos_tmp = cos_wave_shift_times_fi_times_t2{i};
    cos_tmp(1) = cos_tmp(1)/2; cos_tmp(end) = cos_tmp(end)/2;
    int_cos_wave_shift_times_fi_times_t2(i) = sum(cos_tmp)/fs;

    sin_tmp = sin_wave_shift_times_fi_times_t2{i};
    sin_tmp(1) = sin_tmp(1)/2; sin_tmp(end) = sin_tmp(end)/2;
    int_sin_wave_shift_times_fi_times_t2(i) = sum(sin_tmp)/fs;
end

% COMPUTE ENERGY
F_R = sum(int_cos_wave_shift_times_fi);
F_C = -sum(int_sin_wave_shift_times_fi);

E = F_R^2 + F_C^2;

% COMPUTE GRADIENT
grad_F_R = zeros(1,n+1);
grad_F_R(1) = -sum(int_sin_wave_shift_times_fi_times_t)*2*pi;
grad_F_R(2:end) = -int_sin_wave_shift_times_fi(2:end)*2*pi;

grad_F_C = zeros(1,n+1);
grad_F_C(1) = -sum(int_cos_wave_shift_times_fi_times_t)*2*pi;
grad_F_C(2:end) = -int_cos_wave_shift_times_fi(2:end)*2*pi;

grad_E = 2*F_R*grad_F_R + 2*F_C*grad_F_C;
grad_E = grad_E';

% COMPUTE HESSIAN
dwdw_F_R = -sum(int_cos_wave_shift_times_fi_times_t2)*(2*pi)^2;
dwdw_F_C = sum(int_sin_wave_shift_times_fi_times_t2)*(2*pi)^2;

dwdk_F_R = -int_cos_wave_shift_times_fi_times_t(2:end)*(2*pi)^2;
dwdk_F_C = int_sin_wave_shift_times_fi_times_t(2:end)*(2*pi)^2;

dkdk_F_R = -int_cos_wave_shift_times_fi(2:end)*(2*pi)^2;
dkdk_F_C = int_sin_wave_shift_times_fi(2:end)*(2*pi)^2;

hess_E = zeros(n+1,n+1);
for i = 2:n+1
    for j = 2:i
        if i == j
            hess_E(i,j) = 2*grad_F_R(i)^2 + 2*F_R*dkdk_F_R(i-1) ...
                + 2*grad_F_C(i)^2 + 2*F_C*dkdk_F_C(i-1);
        else
            hess_E(i,j) = 2*grad_F_R(i)*grad_F_R(j) + ...
                2*grad_F_C(i)*grad_F_C(j);
        end

        hess_E(j,i) = hess_E(i,j); % symmetric
    end
end

for i = 2:n+1
    hess_E(1,i) = 2*grad_F_R(1)*grad_F_R(i) + 2*F_R*dwdk_F_R(i-1) ...
        + 2*grad_F_C(1)*grad_F_C(i) + 2*F_C*dwdk_F_C(i-1);
    hess_E(i,1) = hess_E(1,i); % symmetric
end

hess_E(1,1) = 2*grad_F_R(1)^2 + 2*F_R*dwdw_F_R + 2*grad_F_C(1)^2 + 2*F_C*dwdw_F_C;

end