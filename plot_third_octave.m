function plot_third_octave(d, e, fs, varargin)
%% plot_third_octave  1/3-octave band power: disturbance vs error
%
%  WHAT IT SHOWS
%  -------------
%  Power (dB) of d and e summed into standard 1/3-octave bands.
%  Each pair of bars (blue = d, orange = e) represents one band.
%  The gap between the two bars is the attenuation your ANC achieves
%  in that frequency region.  A flat, tall blue bar and a short orange
%  bar means strong broadband cancellation.  Isolated tall gaps at a
%  few bands mean purely tonal cancellation — typical for FxLMS on a
%  genset.
%
%  WHAT THE DIFFERENCES MEAN
%  -------------------------
%  Large gap  → ANC is working well at those frequencies.
%  No gap     → ANC is not helping; usually bands between harmonics.
%  e > d      → Controller is adding noise — secondary path mismatch.
%  Gap only at harmonic bands → confirms tonal (not broadband) control.
%
%  INPUTS
%  ------
%  d        : [M × N]  disturbance signal (all mics, or mic 1 row)
%  e        : [M × N]  error signal
%  fs       : scalar   sampling rate (Hz)
%  varargin : optional name-value pairs
%             'mic'     , integer  — which mic row to use (default 1)
%             'f_range' , [flo fhi] — freq limits in Hz (default [20 500])
%             'ss_frac' , scalar  — fraction of signal used as steady
%                                   state window, from the end (default 0.3)
%             'title'   , char    — figure title suffix
%
%  CALL EXAMPLES
%  -------------
%  plot_third_octave(d, e, fs)
%  plot_third_octave(d, e, fs, 'mic', 2, 'f_range', [20 800])
%  plot_third_octave(d, e, fs, 'ss_frac', 0.4, 'title', '1x2x2')

%% ── Parse optional inputs ────────────────────────────────────────────

    % --- NaN-safe: trim to pre-divergence window so spectral
    %     functions never see NaN/Inf after a divergence ---
    [d, e] = nan_safe(d, e);
p = inputParser;
addParameter(p, 'mic',     1,         @isnumeric);
addParameter(p, 'f_range', [20 500],  @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'ss_frac', 0.3,       @(x) isnumeric(x) && x>0 && x<1);
addParameter(p, 'title',   '',        @ischar);
parse(p, varargin{:});
mic     = p.Results.mic;
f_range = p.Results.f_range;
ss_frac = p.Results.ss_frac;
ttl     = p.Results.title;

%% ── Steady-state window ──────────────────────────────────────────────
N      = size(d, 2);
idx_ss = max(1, round(N*(1-ss_frac))) : N;

d_use  = d(mic, idx_ss);
e_use  = e(mic, idx_ss);

%% ── 1/3-octave band power via MATLAB's built-in poctave ─────────────
% poctave (Audio Toolbox) computes fractional-octave band power directly and
% returns the ISO *preferred* centre frequencies (40, 80, 160 ...), which are
% the standard acoustic labels. (Per Prof Gan: prefer validated built-ins.)
%   p  : band power (linear, power in each band)
%   cf : ISO preferred centre frequencies (Hz)
[p_d, cf] = poctave(d_use.', fs, 'BandsPerOctave', 3, ...
                    'FrequencyLimits', f_range);
[p_e, ~ ] = poctave(e_use.', fs, 'BandsPerOctave', 3, ...
                    'FrequencyLimits', f_range);

fc   = cf(:).';                       % row of ISO preferred centres
Pd_b = 10*log10(p_d(:).' + eps);      % band power -> dB
Pe_b = 10*log10(p_e(:).' + eps);
att  = Pd_b - Pe_b;                   % positive = attenuation

%% ── Plot ─────────────────────────────────────────────────────────────
figure('Name', ['1/3-Octave ' ttl], 'Color', 'w');

subplot(2,1,1);
x_pos = 1:numel(fc);
bar(x_pos, [Pd_b; Pe_b].', 'grouped');
set(gca, 'XTick', x_pos, 'XTickLabel', compose('%g', fc), ...
         'XTickLabelRotation', 45, 'FontSize', 9);
ylabel('Power (dB re 1)');
legend('Disturbance d', 'Error e (ANC on)', 'Location', 'northwest');
title(['1/3-Octave Band Power — Mic ' num2str(mic) '  ' ttl]);
grid on;

subplot(2,1,2);
bar(x_pos, att, 'FaceColor', [0.2 0.6 0.3]);
set(gca, 'XTick', x_pos, 'XTickLabel', compose('%g', fc), ...
         'XTickLabelRotation', 45, 'FontSize', 9);
ylabel('Attenuation (dB)');
xlabel('1/3-Octave centre frequency (Hz)');
title('Attenuation = P_d − P_e per band');
yline(0, 'k--');
grid on;

fprintf('[plot_third_octave] Max attenuation: %.1f dB at %.0f Hz\n', ...
        max(att), fc(att == max(att)));
end