function att_h = plot_harmonic_attenuation(d, e, fs, f0, varargin)
%% plot_harmonic_attenuation  Per-harmonic attenuation bar chart
%
%  WHAT IT SHOWS
%  -------------
%  Attenuation (dB) at each integer multiple of the fundamental frequency
%  f0: the generator RPM harmonic (f0, 2f0, 3f0, ...).  Each bar height
%  is how many dB quieter e is than d at that specific harmonic.
%  This is the most direct summary of ANC performance for a tonal source.
%
%  WHAT THE DIFFERENCES MEAN
%  -------------------------
%  Tall bar at f0, shorter at higher harmonics → normal; higher
%      harmonics need more filter taps to model the group delay.
%  Short bar at f0 but tall at 2f0 → possible f0 estimation error;
%      check that f0 is correct.
%  Bar near zero → that harmonic is uncancelled; check coherence
%      between x and d at that frequency (mscohere).
%  Negative bar → controller amplifying that harmonic; secondary
%      path model is wrong in phase at that frequency.
%  Comparing runs: taller bars in run C vs B confirm that a cleaner
%      secondary path estimate directly improves tonal attenuation.
%
%  INPUTS
%  ------
%  d        : [M × N]  disturbance
%  e        : [M × N]  error
%  fs       : scalar   sampling rate (Hz)
%  f0       : scalar   fundamental frequency (Hz) — generator RPM / 60 × poles
%  varargin : optional name-value pairs
%             'mic'      , integer  — mic row (default 1)
%             'n_harm'   , integer  — number of harmonics to plot (default 8)
%             'bw'       , scalar   — half-bandwidth around each harmonic Hz
%                                    (default 5 Hz)
%             'ss_frac'  , scalar   — steady-state fraction from end (default 0.3)
%             'labels'   , cellstr  — extra run labels for grouped bars,
%                                    e.g. {'Run A','Run B','Run C'}
%                                    when supplied, d and e must be cell arrays
%                                    of matching [M×N] matrices
%             'title'    , char     — figure title suffix
%
%  OUTPUT
%  ------
%  att_h    : [1 × n_harm]  attenuation in dB at each harmonic
%
%  CALL EXAMPLES
%  -------------
%  att = plot_harmonic_attenuation(d, e, fs, f0)
%  att = plot_harmonic_attenuation(d, e, fs, f0, 'n_harm', 6, 'bw', 3)
%
%  % Compare multiple runs (grouped bars):
%  att = plot_harmonic_attenuation({d,d,d}, {eA,eB,eC}, fs, f0, ...
%            'labels', {'Run A','Run B','Run C'})

%% ── Parse inputs ─────────────────────────────────────────────────────
p = inputParser;
addParameter(p, 'mic',     1,     @isnumeric);
addParameter(p, 'n_harm',  8,     @isnumeric);
addParameter(p, 'bw',      5,     @isnumeric);
addParameter(p, 'ss_frac', 0.3,   @isnumeric);
addParameter(p, 'labels',  {},    @iscell);
addParameter(p, 'title',   '',    @ischar);
parse(p, varargin{:});
mic     = p.Results.mic;
n_harm  = p.Results.n_harm;
bw      = p.Results.bw;
ss_frac = p.Results.ss_frac;
labels  = p.Results.labels;
ttl     = p.Results.title;

%% ── Handle single vs multi-run input ─────────────────────────────────
if iscell(d)
    n_runs = numel(d);
else
    d      = {d};
    e      = {e};
    n_runs = 1;
    if isempty(labels), labels = {'ANC'}; end
end
if isempty(labels)
    labels = arrayfun(@(k) sprintf('Run %d',k), 1:n_runs, 'uni', 0);
end

%% ── Harmonics ────────────────────────────────────────────────────────
harm = f0 * (1:n_harm);

%% ── Compute attenuation per run ──────────────────────────────────────
att_all = zeros(n_runs, n_harm);

for r = 1:n_runs
    N_r    = size(d{r}, 2);
    idx_ss = max(1, round(N_r*(1-ss_frac))) : N_r;
    d_use  = d{r}(mic, idx_ss);
    e_use  = e{r}(mic, idx_ss);

    nfft   = min(8192, numel(d_use));
    [Pxx_d, F] = pwelch(d_use, hann(nfft), floor(nfft/2), nfft, fs);
    [Pxx_e, ~] = pwelch(e_use, hann(nfft), floor(nfft/2), nfft, fs);

    for h = 1:n_harm
        mask = abs(F - harm(h)) <= bw;
        if ~any(mask)
            att_all(r,h) = NaN;
        else
            Pd = mean(Pxx_d(mask));
            Pe = mean(Pxx_e(mask));
            att_all(r,h) = 10*log10((Pd + eps) / (Pe + eps));
        end
    end
end

att_h = att_all(1,:);   % return first run for single-run case

%% ── Compute harmonic-weighted average attenuation ────────────────────
% Weight each harmonic by its share of total disturbance power
for r = 1:n_runs
    N_r    = size(d{r}, 2);
    idx_ss = max(1, round(N_r*(1-ss_frac))) : N_r;
    d_use  = d{r}(mic, idx_ss);
    nfft   = min(8192, numel(d_use));
    [Pxx_d, F] = pwelch(d_use, hann(nfft), floor(nfft/2), nfft, fs);
    harm_pow = zeros(1,n_harm);
    for h = 1:n_harm
        mask = abs(F - harm(h)) <= bw;
        if any(mask), harm_pow(h) = mean(Pxx_d(mask)); end
    end
    w_avg = sum(harm_pow .* att_all(r,:), 'omitnan') / (sum(harm_pow) + eps);
    fprintf('[plot_harmonic_attenuation] %s — weighted avg attenuation: %.1f dB\n', ...
            labels{r}, w_avg);
end

%% ── Plot ─────────────────────────────────────────────────────────────
figure('Name', ['Harmonic Attenuation ' ttl], 'Color', 'w');

if n_runs == 1
    bar(1:n_harm, att_h, 'FaceColor', [0.2 0.5 0.8]);
else
    bar(1:n_harm, att_all.');   % grouped bars, one group per harmonic
    legend(labels, 'Location', 'northeast');
end

set(gca, 'XTick', 1:n_harm, ...
         'XTickLabel', arrayfun(@(h) sprintf('%.0f Hz', h), harm, 'uni', 0));
xlabel('Harmonic');
ylabel('Attenuation (dB)');
title(['Per-Harmonic Attenuation   f_0 = ' num2str(f0,'%.1f') ' Hz   ' ttl]);
yline(0, 'k--', 'LineWidth', 1);
grid on;
end