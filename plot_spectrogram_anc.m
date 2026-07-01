function plot_spectrogram_anc(d, e, fs, varargin)
%% plot_spectrogram  Side-by-side spectrogram: disturbance vs error
%
%  WHAT IT SHOWS
%  -------------
%  A time-frequency map of d (left) and e (right).  Colour = power in dB.
%  The x-axis is time, the y-axis is frequency, so you can see BOTH:
%    - which frequencies carry energy (horizontal bright stripes = harmonics)
%    - when those frequencies disappear (stripes fading = ANC converging)
%  This is the only plot that shows convergence dynamics in time.
%
%  WHAT THE DIFFERENCES MEAN
%  -------------------------
%  Bright stripe in d, faded in e   → ANC cancelled that harmonic.
%  Stripe fades slowly (>3 s)       → slow convergence; mu may be too small
%                                     or Lw too large.
%  Stripe fades then reappears      → filter diverging or path changed.
%  New stripes appear in e not in d → nonlinear distortion or aliasing;
%                                     check loudspeaker clipping.
%  Both panels identical            → ANC has no effect; check wiring or
%                                     that e is not just a copy of d.
%
%  INPUTS
%  ------
%  d        : [M × N]  disturbance
%  e        : [M × N]  error
%  fs       : scalar   sampling rate (Hz)
%  varargin : optional name-value pairs
%             'mic'    , integer     — mic row to plot (default 1)
%             'f_max'  , scalar      — max frequency to display Hz (default 500)
%             'clim'   , [lo hi]     — colour axis limits dB (default [-120 -40])
%             'win'    , integer     — STFT window length in samples (default 2048)
%             'title'  , char        — figure title suffix
%
%  CALL EXAMPLES
%  -------------
%  plot_spectrogram(d, e, fs)
%  plot_spectrogram(d, e, fs, 'mic', 2, 'f_max', 800)
%  plot_spectrogram(d, e, fs, 'clim', [-100 -30], 'win', 1024)

%% ── Parse inputs ─────────────────────────────────────────────────────
p = inputParser;
addParameter(p, 'mic',   1,           @isnumeric);
addParameter(p, 'f_max', 500,         @isnumeric);
addParameter(p, 'clim',  [-120 -40],  @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'win',   2048,        @isnumeric);
addParameter(p, 'title', '',          @ischar);
parse(p, varargin{:});
mic   = p.Results.mic;
f_max = p.Results.f_max;
clim_ = p.Results.clim;
Nwin  = p.Results.win;
ttl   = p.Results.title;

hop  = floor(Nwin / 2);
nfft = Nwin;

d_use = d(mic, :);
e_use = e(mic, :);

%% ── Plot ─────────────────────────────────────────────────────────────
figure('Name', ['Spectrogram ' ttl], 'Color', 'w');

subplot(1,2,1);
spectrogram(d_use, hann(Nwin), hop, nfft, fs, 'yaxis');
title(['Disturbance  d  (mic ' num2str(mic) ')']);
ylim([0 f_max/1000]);   % spectrogram yaxis is in kHz
clim(clim_);
colormap(gca, 'hot');
xlabel('Time (s)');  ylabel('Frequency (Hz)');

subplot(1,2,2);
spectrogram(e_use, hann(Nwin), hop, nfft, fs, 'yaxis');
title(['Error  e  — ANC on  (mic ' num2str(mic) ')   ' ttl]);
ylim([0 f_max/1000]);
clim(clim_);
colormap(gca, 'hot');
xlabel('Time (s)');  ylabel('Frequency (Hz)');

% Shared colourbar on right
cb = colorbar('peer', gca);
ylabel(cb, 'Power (dB)');

sgtitle(['Spectrogram before/after ANC   ' ttl], 'FontWeight', 'bold');

%% ── Print rough convergence times (stripes crossing -60 dB) ─────────
% Compute power vs time at each 1-second window to find when it drops
N      = size(d,2);
t_full = (0:N-1)/fs;
win_s  = round(fs);   % 1-second averaging window
steps  = floor(N / win_s);
t_win  = (0.5 : steps-0.5) * win_s / fs;

pwr_d = zeros(1, steps);
pwr_e = zeros(1, steps);
for s = 1:steps
    idx = (s-1)*win_s+1 : s*win_s;
    pwr_d(s) = 10*log10(mean(d_use(idx).^2) + eps);
    pwr_e(s) = 10*log10(mean(e_use(idx).^2) + eps);
end
target_dB = pwr_d(1) - 10;   % 10 dB below initial disturbance
idx_conv  = find(pwr_e < target_dB, 1);
if ~isempty(idx_conv)
    fprintf('[plot_spectrogram] Approx convergence (−10 dB) at t ≈ %.1f s\n', ...
            t_win(idx_conv));
else
    fprintf('[plot_spectrogram] Signal did not drop 10 dB below initial level.\n');
end
end