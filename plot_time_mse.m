function plot_time_mse(cfg, e, n_div, varargin)
%% plot_time_mse  Time-domain MSE convergence curve
%
%  WHAT IT SHOWS
%  -------------
%  Mean squared error (averaged across all error mics) in dB, smoothed
%  with a moving average, plotted against time. Shows how quickly the
%  adaptive filter converges and whether/where it diverges.
%
%  INPUTS
%  ------
%  cfg      : struct — needs fs, T, K, M, muw
%  e        : [M x N]  error signal (output of MCANC_FxLMS2)
%  n_div    : scalar — divergence sample index, NaN if stable
%  varargin : optional name-value pairs
%             'label' , char — run label for titles/legend (default 'ANC Run')
%
%  CALL EXAMPLE
%  ------------
%  plot_time_mse(cfg, e, n_div, 'label', 'Run C');

p = inputParser;
addParameter(p, 'label', 'ANC Run', @ischar);
parse(p, varargin{:});
lbl = p.Results.label;

N = cfg.fs * cfg.T;
t = (0:N-1) / cfg.fs;

figure('Name', [lbl ' — Fig 1: MSE Convergence'], 'Color', 'w');
hold on; grid on;

mse_mean = mean(e.^2, 1, 'omitnan');
mse_db   = 10*log10(movmean(mse_mean, round(0.1*cfg.fs), 'omitnan') + eps);
plot(t, mse_db, 'Color', [0.18 0.45 0.85], 'LineWidth', 1.6, 'DisplayName', lbl);

if ~isnan(n_div)
    xline(n_div/cfg.fs, 'r:', 'LineWidth', 1.5, 'DisplayName', 'Divergence');
end

xlabel('Time (s)');
ylabel('Mean error power (dB)');
title(sprintf('%s — MSE Convergence  (K=%d, M=%d, \\mu=%.1e)', lbl, cfg.K, cfg.M, cfg.muw));
legend('Location', 'best');
xlim([0, min(10, cfg.T)]);

end