clear
close all
clc

rng(0)   
example = 2;

%%
if example == 1
    %% EXAMPLE 1: 1 SEGMENT OF DATA
    load('1segment_example.mat') 

    N = length(S);

    % add artificial gap of 0
    S_cell = {}; S_cell{1} = S(1:floor(N/2)); S_cell{2} = S(floor(N/2)+1:end); 
    samp_shift = 0;
elseif example == 2
    %% EXAMPLE 2: MULTIPLE SEGMENTS OF DATA
    load('multsegment_example.mat') 
    
    % convert S (gap padded with NaNs) to cell array
    S_cell = convert_vector_to_cellarray(S, N, samp_shift); 
end

%% RUN ALGORITHM 1
K = 10;     % number of harmonics to fit
w0 = 150.6; % initial guess  
tic
[w, d, ~, t] = newton_rand_init(w0, 5, 25, 5000, 1000, S_cell, fs, 1e-8);
[B_est, A_est, alp, t_vec] = remove_artifact(S_cell, t, fs, K, w, d); % if only using algorithm 2, skip this step
toc
B_est_vec = convert_cellarray_to_vector(B_est, samp_shift, 0); % convert to vectors
A_est_vec = convert_cellarray_to_vector(A_est, zeros(length(S)-1), nan);

%% RUN ALGORITHM 2 
tic
% initialize with outputs of newton_rand_init
[w_refine, d_refine, ~, B_est_refine, A_est_refine, alp_refine, t_vec_refine] = ...
    newton_refinement_using_g(w, d, 1000, S_cell, t, fs, K, 1e-8); 
toc
B_est_vec_refine = convert_cellarray_to_vector(B_est_refine, samp_shift, 0); % convert to vectors
A_est_vec_refine = convert_cellarray_to_vector(A_est_refine, zeros(length(S)-1), nan);

%% PLOT RESULTS
% RECONSTRUCTED ARTIFACT, MODDED TIME PLOTS
figure
S = convert_cellarray_to_vector(S_cell, zeros(1,length(S_cell)-1), 0);
subplot(1,2,1)
plot(mod(t_vec_true,1/freq_true), S, '.')
hold on
plot(mod(t_vec,1/w), A_est_vec, '.')
title('Algorithm 1')
legend('observed signal', 'reconstructed artifact')
xlabel('time modded by period')

subplot(1,2,2)
plot(mod(t_vec_true,1/freq_true), S, '.')
hold on
plot(mod(t_vec_refine,1/w_refine), A_est_vec_refine, '.')
title('Algorithm 2')
legend('observed signal', 'reconstructed artifact')
xlabel('time modded by period')

%% ERROR IN TIME DOMAIN
figure
tiledlayout(1,2)
ax1 = nexttile;
t_plot = (0:length(B)-1)/fs;
plot(t_plot, B_est_vec - B)
title('Error in Signal Recovered by Algorithm 1')
xlabel('time')

ax2 = nexttile;
plot(t_plot, B_est_vec_refine - B)
title('Error in Signal Recovered by Algorithm 2')
xlabel('time')

linkaxes([ax1, ax2], 'xy')

%% POWER SPECTRUM USING SPECTROGRAM
S = convert_cellarray_to_vector(S_cell, samp_shift, 0);
figure
subplot(1,4,1)
spectrogram(S, [], [], [], fs, 'yaxis')
title('observed signal')

subplot(1,4,2)
spectrogram(B, [], [], [], fs, 'yaxis')
title('true underlying signal')

subplot(1,4,3)
spectrogram(B_est_vec, [], [], [], fs, 'yaxis')
title('algorithm 1')

subplot(1,4,4)
spectrogram(B_est_vec_refine, [], [], [], fs, 'yaxis')
title('algorithm 2')

%% POWER SPECTRUM USING PWELCH
figure
tiledlayout(1,4)
ax1 = nexttile;
pwelch(S,[],[],0:fs/2,fs)
title('observed signal')

ax2 = nexttile;
pwelch(B,[],[],0:fs/2,fs)
title('true underlying signal')

ax3 = nexttile;
pwelch(B_est_vec,[],[],0:fs/2,fs)
title('algorithm 1')

ax4 = nexttile;
pwelch(B_est_vec_refine,[],[],0:fs/2,fs)
title('algorithm 2')

linkaxes([ax1 ax2 ax3 ax4], 'xy')