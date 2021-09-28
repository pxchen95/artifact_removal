function [A, t, t_shift, amp] = create_artifact(ampC, ampS, numHarmonics, freq, N, samp_shift, fs)
% Create samples of periodic artifact of the form: 
% A(t) = amp_0 + sum_i^K amp_{i,c}*cos(2*pi*w*i*t) + amp_{i,s}*sin(2*pi*w*i*t))
% with numSegments-1 gaps and K = numHarmonics 
%
% INPUTS:
%     ampC:         scalar, amplitude of cos(2*pi*t*freq)
%     ampS:         scalar, amplitude of sin(2*pi*t*freq)
%     numHarmonics: scalar, number of harmonics
%     freq:         scalar, fundamental frequency
%     N:            1 x numSegments vector, N(i) = length of segment i in samples 
%     samp_shift:   1 x numSegments vector, samp_shift(i) = length of gap b/t segments i,i+1 in samples
%     fs:           scalar, sample rate
%
% OUTPUTS: 
%     A:       1 x numSegments cell array, A{i} = 1 x N_i, artifact samples of segment i
%     t:       1 x numSegments cell array, t{i} = 1 x N_i, unshifted times that segment i is sampled at, in [0, T_i]
%     t_shift: 1 x numSegments cell array, t{i} = 1 x N_i, shifted times that segment i is sampled at
%     amp:     1 x 2*numHarmonics, amp(1:numHarmonics) = amplitude of cosine waves, amp(numHarmonics+1:end) = amplitude of sine waves

% COMPUTE SAMPLE TIMES
t = {}; t_shift = {};
numSegments = length(N);
sumT = 0;
sum_shift = 0; 
for i = 1:numSegments
    N_i = N(i);    % number of samples in segment i
    T_i = N_i/fs;  % length of segment i in time

    tshift = sumT + sum_shift; % segment i starts sampling at time tshift

    t_i = linspace(0, T_i, N_i+1);  % sample times for segment i
    t_i(end) = []; 
    t{i} = t_i;         % unshifted sample times, \in [0, T_i]
    t_i = t_i + tshift; % shifted sample times, \in [t_shift, t_shift + t_i]
    t_shift{i} = t_i;

    if i < numSegments
        samp_shift_i = samp_shift(i); % number of samples segment i is shifted over

        sumT = sumT + T_i; % update sums
        sum_shift = sum_shift + samp_shift_i/fs;
    end
end

% CREATE ARTIFACT
amp = zeros(2*numHarmonics, 1);
A = {}; 
for k = 1:numHarmonics % sum harmonics
    if k == 1
        amp(k) = ampC;
        amp(numHarmonics + k) = ampS;
    else % amplitudes of harmonics decay + some randomness
        amp(k) = ampC * rand(1) * 2 * 2^(-(k-1));                % cosine amplitudes
        amp(numHarmonics + k) = ampS * rand(1) * 2 * 2^(-(k-1)); % sine amplitudes
    end
    A_func = @(t) amp(numHarmonics + k)*sin(2*pi*t*k*freq) + amp(k)*cos(2*pi*t*k*freq); % kth harmonic
   
    for i = 1:numSegments
        t_i = t_shift{i};
        if k == 1
            A{i} = A_func(t_i);
        else
            A{i} = A{i} + A_func(t_i); % artifact samples in segment i
        end
    end
end

end