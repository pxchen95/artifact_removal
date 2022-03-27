clear
close all
clc

%% CHOOSE EXAMPLE TO RUN
example = 3;

if example == 1 % single segment, simulated artifact, no underlying signal
    load('paper_examples/example1_SingleSegmentArtifactOnly.mat')
    A_vec_0   = A_true{1}; % vectors formatted with appropriate padding for plotting/error calcs
    A_vec_nan = A_vec_0;
    B_vec_nan = B_true{1};
    S_vec_0   = S{1}; 
    S_vec_nan = S_vec_0;
elseif example == 2 % single segment, simulated artifact with a chirp
    load('paper_examples/example2_SingleSegmentChirp.mat')
    A_vec_0   = A_true{1}; % vectors formatted with appropriate padding for plotting/error calcs
    A_vec_nan = A_vec_0;
    B_vec_nan = B_true{1};
    B_vec_0   = B_vec_nan;
    S_vec_0   = S{1}; 
    S_vec_nan = S_vec_0;
elseif example == 3 % multiple segments, aliased simulated artifact with simulated underlying signal
    load('paper_examples/example3_ManySegmentsAliased.mat')
    A_vec_0   = convert_cellarray_to_vector(A_true, samp_shift, 0); % vectors formatted with appropriate padding for plotting/error calcs
    A_vec_nan = convert_cellarray_to_vector(A_true, samp_shift, nan);
    B_vec_0   = convert_cellarray_to_vector(B_true, samp_shift, 0);
    B_vec_nan = convert_cellarray_to_vector(B_true, samp_shift, nan);
    S_vec_0   = convert_cellarray_to_vector(S, samp_shift, 0);
    S_vec_nan = convert_cellarray_to_vector(S, samp_shift, nan);
end

%%
rng(0)
K = 5; % # of harmonics to fit

%% RUN INITIALIZATION ALGORITHM
tic
[w_NA, d_NA, ~, t] = newton_rand_init(150.6, 5, 25, 5000, 1000, S, fs, 1e-8);
toc

%% RUN ALGORITHM 1
tic
[w_est, d_est, ~, B_est, A_est, ~, t_vec] = ...
    newton_refinement_using_g(w_NA, d_NA, 1000, S, t, fs, K, 1e-8); 
toc

%% CONVERT CELL ARRAYS TO VECTORS 
B_est_vec_0 = convert_cellarray_to_vector(B_est, samp_shift, 0);     % padded with 0s
B_est_vec_nan = convert_cellarray_to_vector(B_est, samp_shift, nan); % padded with nans
A_est_vec = convert_cellarray_to_vector(A_est, zeros(length(S)-1), nan);

%% DISPLAY ERRORS
ind = ~isnan(B_est_vec_nan); % indices of non-NaN entries
disp(['Relative Error of Frequency Estimate: ' ...
    num2str(abs(w_est - freq_true)/freq_true*100) '%'])
disp(['Relative RMSE of Reconstructed Artifact: ' ...
    num2str(norm(A_est_vec - A_vec_nan(ind))/norm(A_vec_nan(ind))*100) '%'])
if example == 1
    disp(['RMSE of Recovered Signal: ' num2str(norm(B_est_vec_0))])
else
    disp(['Relative RMSE of Recovered Signal: ' ...
        num2str(norm(B_est_vec_0(ind) - B_vec_nan(ind))/norm(B_vec_nan(ind))*100) '%'])
end

%% PLOT RESULTS
if example == 1
    %% PLOT ENERGY
    freq_center = freq_true; freq_center_zoom = w_NA;
    freq_width = 0.5; freq_width_zoom = 5e-6;
    numPoints = 1000;
    OmegaE = linspace(freq_center-freq_width, freq_center+freq_width, numPoints);
    OmegaE_zoom = linspace(freq_center_zoom-freq_width_zoom, freq_center_zoom+freq_width_zoom, numPoints);
    E = zeros(1,numPoints); E_zoom = E;
    for i = 1:numPoints
        E(i) = calc_E_multtshifts(OmegaE(i), [], S, {t_vec_true}, fs);
        E_zoom(i) = calc_E_multtshifts(OmegaE_zoom(i), [], S, {t_vec_true}, fs);
    end

    E_freq_true = calc_E_multtshifts(freq_true, [], S, {t_vec_true}, fs);
    E_w_NA = calc_E_multtshifts(w_NA, [], S, {t_vec_true}, fs);

    %%
    figure
    subplot(1,2,1)
    plot(OmegaE, E, 'DisplayName', 'E(\omega)')
    hold on
    plot(freq_true, E_freq_true, '*k', 'DisplayName', 'True Frequency')
    plot(w_NA, E_w_NA, 'or', 'MarkerFaceColor', 'r', 'DisplayName', 'Algorithm 2 Estimate')
    legend('AutoUpdate', 'off', 'Location', 'SouthWest')
    plot(freq_true, E_freq_true, '*k')
    xlim([OmegaE(1), OmegaE(end)])
    xlabel('$\omega$', 'interpreter', 'latex')
    ylabel('E')
    grid on
    
    subplot(1,2,2)
    plot(OmegaE_zoom, E_zoom, 'DisplayName', 'E(\omega)')
    hold on
    plot(freq_true, E_freq_true, '*k', 'DisplayName', 'True Frequency')
    plot(w_NA, E_w_NA, 'or', 'MarkerFaceColor', 'r', 'DisplayName', 'Algorithm 2 Estimate')
    legend('AutoUpdate', 'off', 'Location', 'SouthWest')
    xlim([OmegaE_zoom(1), OmegaE_zoom(end)])
    xlabel('$\omega$', 'interpreter', 'latex')
    ylabel('E')
    grid on
    
    %% PLOT LEAST SQUARES OBJECTIVE FUNCTION (g)
    numPoints = 100;
    OmegaG = linspace(freq_center-freq_width, freq_center+freq_width, numPoints);
    OmegaG_zoom = linspace(freq_center-freq_width_zoom, freq_center+freq_width_zoom, numPoints);
    g = zeros(size(OmegaG)); g_zoom = g;
    for i = 1:numPoints
        [~, ~, ~, ~, ~, g(i)] = ...
            remove_artifact_ver_g(S, {t_vec_true}, fs, 15, OmegaG(i), []);
        [~, ~, ~, ~, ~, g_zoom(i)] = ...
            remove_artifact_ver_g(S, {t_vec_true}, fs, 15, OmegaG_zoom(i), []);
    end

    %%
    [~, ~, ~, ~, ~, g_freq_true] = ...
                remove_artifact_ver_g(S, {t_vec_true}, fs, 15, freq_true, []);
    [~, ~, ~, ~, ~, g_w_NA] = ...
                remove_artifact_ver_g(S, {t_vec_true}, fs, 15, w_NA, []);       
    [~, ~, ~, ~, ~, g_w_refine] = ...
                remove_artifact_ver_g(S, {t_vec_true}, fs, 15, w_est, []);      

    %%
    figure
    subplot(1,2,1)
    plot(OmegaG, g, 'DisplayName', 'g(\omega)')
    hold on
    plot(freq_true, g_freq_true, '*k', 'DisplayName', 'True Frequency')
    plot(w_est, g_w_refine, 'gs', 'MarkerFaceColor', 'g', 'DisplayName', 'Algorithm 1 Estimate')
    plot(w_NA, g_w_NA, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Algorithm 2 Estimate')
    xlim([OmegaG(1) OmegaG(end)])
    legend('AutoUpdate', 'off', 'Location', 'SouthWest')
    plot(freq_true, g_freq_true, '*k')
    grid on
    xlabel('$\omega$', 'interpreter', 'latex')
    ylabel('g')
    
    subplot(1,2,2)
    plot(OmegaG_zoom, g_zoom, 'DisplayName', 'g(\omega)')
    hold on
    plot(freq_true, g_freq_true, '*k', 'DisplayName', 'True Frequency')
    plot(w_est, g_w_refine, 'gs', 'MarkerFaceColor', 'g', 'DisplayName', 'Algorithm 1 Estimate')
    plot(w_NA, g_w_NA, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Algorithm 2 Estimate')
    xlim([OmegaG_zoom(1) OmegaG_zoom(end)])
    legend('AutoUpdate', 'off', 'Location', 'SouthWest')
    plot(freq_true, g_freq_true, '*k')
    grid on
    xlabel('$\omega$', 'interpreter', 'latex')
    ylabel('g')
end

%% PLOT RECONSTRUCTED ARTIFACT
figure
hold on
if example ~= 1, plot(mod(t_vec_true,1/w_est), S_vec_nan, 's', 'DisplayName', 'Observed Signal'); end
plot(mod(t_vec_true,1/freq_true), A_vec_nan, '+', 'DisplayName', 'True Artifact')
plot(mod(t_vec,1/w_est), A_est_vec, '.', 'DisplayName', 'Reconstructed Artifact')
xlim([0 max(mod(t_vec,1/w_est))])
legend()
xlabel('Time Modulo the Period (s)')
ylabel('Voltage')

%% PLOT RECOVERED SIGNALS
figure
hold on
if example == 3, plot(S_vec_nan, 'Color', [0.9290 0.6940 0.1250], 'DisplayName', 'Observed Signal'); end
plot(B_vec_nan, 'Color', [0 0.4470 0.7410], 'DisplayName', 'True Underlying Signal')
plot(B_est_vec_nan, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'Recovered Signal')
legend('AutoUpdate', 'off', 'Location', 'SouthEast')
if example == 1, plot(B_vec_nan, 'Color', [0, 0.4470, 0.7410]); end
ylabel('Voltage')
xlabel('Sample Number')
xlim([1, length(B_vec_nan)])

%% PWELCH
if example == 1 || example == 2
    figure
    if example == 1, numTiles = 2; else, numTiles = 3; end
    tiledlayout(1,numTiles)
    ax1 = nexttile;
    pwelch(S_vec_0,[],[],0:fs/2,fs)
    title('Observed Signal')

    ax2 = nexttile;
    pwelch(B_est_vec_0,[],[],0:fs/2,fs)
    if example == 1, ylim([-280,0]); end
    title('Recovered Signal')
    ylabel('')

    if example ~= 1
        ax3 = nexttile;
        pwelch(B_vec_0,[],[],0:fs/2,fs)
        title('True Underlying Signal')
        ylabel('')
        linkaxes([ax1, ax2, ax3], 'xy')
    end
end

%% SPECTROGRAM
figure
if example == 1
    colorbar_limits = [-50 0];
    window = 500;
    noverlap = 125;
    nfft = 500;
elseif example == 2
    colorbar_limits = [-65 0];
    window = 325;
    noverlap = 250;
    nfft = 500;
elseif example == 3
    colorbar_limits = [-40 0];
    window = 125;
    noverlap = 25;
    nfft = 500;
end
if example == 1, numPlots = 2; else, numPlots = 3; end

subplot(1,numPlots,1)
spectrogram(S_vec_0, window, noverlap, nfft, fs, 'yaxis')
title('Observed Signal')
caxis(colorbar_limits)

subplot(1,numPlots,2)
spectrogram(B_est_vec_0, window, noverlap, nfft, fs, 'yaxis')
title('Recovered Signal')
caxis(colorbar_limits)

if example ~= 1
    subplot(1,numPlots,3)
    spectrogram(B_vec_0, window, noverlap, nfft, fs, 'yaxis')
    title('True Underlying Signal')
    caxis(colorbar_limits)
end