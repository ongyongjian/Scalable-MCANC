function plot_nr_vs_achieved(cfg, d, e, idx_ss, coherence_results, varargin)
%% plot_nr_vs_achieved  Achieved attenuation vs coherence-derived NR bound
%
%  WHAT IT SHOWS
%  -------------
%  Overlays the theoretical maximum attenuation (the coherence-based NR
%  ceiling from check_coherence_NR.m) against what MCANC_FxLMS2 actually
%  achieved, per error mic. The gap between the two curves shows how
%  much headroom remains — whether from filter length, step size, or
%  secondary-path identification quality.
%
%  INPUTS
%  ------
%  cfg               : struct — needs fs, M, f_hi
%  d                 : [M x N]  disturbance
%  e                 : [M x N]  error signal
%  idx_ss            : steady-state sample index range to analyse
%  coherence_results : struct array from check_coherence_NR.m
%  varargin          : optional name-value pairs
%                       'label'      , char    — run label
%                       'harm_freqs' , vector  — harmonic markers (Hz)
%
%  CALL EXAMPLE
%  ------------
%  plot_nr_vs_achieved(cfg, d, e, idx_ss, coherence_results, ...
%                       'label', 'Run C', 'harm_freqs', harm_freqs);

p = inputParser;
addParameter(p, 'label',      'ANC Run', @ischar);
addParameter(p, 'harm_freqs', [],        @isnumeric);
parse(p, varargin{:});
lbl        = p.Results.label;
harm_freqs = p.Results.harm_freqs;

if isempty(coherence_results)
    fprintf('[plot_nr_vs_achieved] No coherence_results supplied — skipping.\n');
    return;
end

win_psd  = hamming(4096);
nov_psd  = 2048;
nfft_psd = 8192;

for m = 1:cfg.M
    [PSD_d, f_p] = pwelch(d(m, idx_ss), win_psd, nov_psd, nfft_psd, cfg.fs);
    [PSD_e, ~]   = pwelch(e(m, idx_ss), win_psd, nov_psd, nfft_psd, cfg.fs);
    atten_achieved = 10*log10(PSD_d ./ (PSD_e + eps));

    figure('Name', [lbl sprintf(' — Fig 6: NR Bound vs Achieved Mic %d', m)], 'Color', 'w');
    hold on; grid on;

    plot(coherence_results(m).f, coherence_results(m).NR_dB, 'k--', ...
         'LineWidth', 1.6, 'DisplayName', 'Coherence NR bound (ceiling)');
    plot(f_p, atten_achieved, 'Color', [0.10 0.62 0.18], 'LineWidth', 1.4, ...
         'DisplayName', [lbl ' (achieved)']);

    for h = 1:numel(harm_freqs)
        xline(harm_freqs(h), 'k:', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    end

    xlabel('Frequency (Hz)'); ylabel('Attenuation (dB)');
    title(sprintf('%s Mic %d — Achieved vs Theoretical Ceiling', lbl, m));
    legend('Location', 'best');
    xlim([0, cfg.f_hi]);

    fprintf('[plot_nr_vs_achieved] Mic %d: median achieved = %.1f dB | median ceiling = %.1f dB\n', ...
            m, median(atten_achieved(isfinite(atten_achieved))), coherence_results(m).NR_median_dB);
end

end