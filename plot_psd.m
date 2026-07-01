function plot_psd(cfg, d, e, idx_ss, varargin)
%% plot_psd  PSD before/after ANC, per error microphone
%
%  Power spectral density of d (before ANC) overlaid with e (after ANC),
%  in the steady-state window, for each error mic. Dotted vertical lines
%  mark the detected harmonics.
%
%  INPUTS
%  ------
%  cfg      : struct — needs fs, M, f_hi
%  d        : [M x N]  disturbance
%  e        : [M x N]  error signal
%  idx_ss   : steady-state sample index range to analyse
%  varargin : optional name-value pairs
%             'label'      , char    — run label for titles/legend
%             'harm_freqs' , vector  — harmonic frequencies to mark (Hz)
%
%  CALL EXAMPLE
%  ------------
%  plot_psd(cfg, d, e, idx_ss, 'label', 'Run C', 'harm_freqs', harm_freqs);

p = inputParser;
addParameter(p, 'label',      'ANC Run', @ischar);
addParameter(p, 'harm_freqs', [],        @isnumeric);
parse(p, varargin{:});
lbl        = p.Results.label;
harm_freqs = p.Results.harm_freqs;

win_psd  = hamming(4096);
nov_psd  = 2048;
nfft_psd = 8192;

for m = 1:cfg.M
    figure('Name', [lbl sprintf(' — Fig 2: PSD Mic %d', m)], 'Color', 'w');
    hold on; grid on;

    [PSD_d, f_p] = pwelch(d(m, idx_ss), win_psd, nov_psd, nfft_psd, cfg.fs);
    [PSD_e, ~]   = pwelch(e(m, idx_ss), win_psd, nov_psd, nfft_psd, cfg.fs);

    atten_m = 20*log10(rms(d(m,idx_ss)) / (rms(e(m,idx_ss)) + eps));

    plot(f_p, 10*log10(PSD_d + eps), 'k', 'LineWidth', 2.2, 'DisplayName', 'Before ANC (d)');
    plot(f_p, 10*log10(PSD_e + eps), 'Color', [0.85 0.25 0.15], 'LineWidth', 1.5, ...
         'DisplayName', sprintf('%s (%+.1f dB)', lbl, atten_m));

    for h = 1:numel(harm_freqs)
        xline(harm_freqs(h), 'k:', 'LineWidth', 0.6, 'HandleVisibility', 'off');
    end

    xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
    title(sprintf('%s — Mic %d Steady-State PSD', lbl, m));
    legend('Location', 'best');
    xlim([0, cfg.f_hi]);
end

end