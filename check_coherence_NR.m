function results = check_coherence_NR(x, d, fs, varargin)
%% check_coherence_NR  Coherence-derived noise-reduction (NR) bound
%                      General J x K x M.
%
%  WHAT THIS ANSWERS
%  ------------------
%  Before running any ANC: what is the theoretical maximum noise reduction
%  at each frequency, given the reference signal(s) and the disturbance?
%  Per Kuo & Morgan, "Active Noise Control Systems", Section 3.1.2.
%
%    Single reference (J = 1), per error mic m:
%      C_dx(w)  = |S_dx(w)|^2 / ( S_dd(w) * S_xx(w) )         Eq.(3.1.4)
%      NR(dB)   = -10*log10( 1 - C_dx(w) )                     p.57
%
%    Multiple INDEPENDENT references (J > 1), per error mic m:
%      Uses the MULTIPLE COHERENCE function
%        C_{d|x}(w) = p^H(w) * Sxx^{-1}(w) * p(w) / S_dd(w)
%      where
%        p(w)   = [S_{d,x1}(w); ...; S_{d,xJ}(w)]   (J x 1 cross-spectra)
%        Sxx(w) = J x J reference cross-spectral matrix
%      This collapses to the ordinary C_dx above when J = 1, so a single-
%      reference run is unchanged. NR(dB) = -10*log10(1 - C_{d|x}).
%
%  Each error mic m = 1..M is treated independently (its own disturbance
%  d(m,:) against the shared set of J references).
%
%  INPUTS
%  ------
%  x        : [J x N]   reference signal(s), one row per reference sensor
%  d        : [M x N]   disturbance at each of the M error microphones
%  fs       : scalar    sampling rate (Hz)
%  varargin : optional name-value pairs
%             'nfft'    , integer  — FFT length (default 8192)
%             'window'  , vector   — window function (default hann(4096))
%             'noverlap', integer  — overlap samples (default nfft/2)
%             'f_range' , [flo fhi] — band to report/plot (default [0 fs/2])
%             'oct_band', logical  — also compute 1/3-octave summary (default true)
%             'plot'    , logical  — show figures (default true)
%             'title'   , char     — figure title suffix
%
%  OUTPUT  (struct array, one entry per mic)
%  -------
%  results(m).f             : [Nf x 1]  frequency axis (Hz)
%  results(m).C_dx          : [Nf x 1]  (multiple) coherence at mic m
%  results(m).NR_dB         : [Nf x 1]  max achievable NR (dB) at mic m
%  results(m).fc_oct        : [1 x Nb]  1/3-oct centre frequencies
%  results(m).C_dx_oct      : [1 x Nb]  1/3-oct band-averaged coherence
%  results(m).NR_oct_dB     : [1 x Nb]  1/3-oct band NR bound (dB)
%  results(m).NR_median_dB  : scalar    median NR bound across the band
%  results(m).NR_best_dB    : scalar    best (highest) NR bound, single bin
%
%  CALL EXAMPLE
%  ------------
%  results = check_coherence_NR(x, d, cfg.fs, 'f_range', [cfg.f_lo cfg.f_hi]);

%% ── Parse inputs ─────────────────────────────────────────────────────
p = inputParser;
addParameter(p, 'nfft',     8192,        @isnumeric);
addParameter(p, 'window',   hann(4096),  @isnumeric);
addParameter(p, 'noverlap', 2048,        @isnumeric);
addParameter(p, 'f_range',  [0 fs/2],    @(v) isnumeric(v) && numel(v)==2);
addParameter(p, 'oct_band', true,        @islogical);
addParameter(p, 'plot',     true,        @islogical);
addParameter(p, 'title',    '',          @ischar);
parse(p, varargin{:});
nfft     = p.Results.nfft;
win      = p.Results.window;
noverlap = p.Results.noverlap;
f_range  = p.Results.f_range;
do_oct   = p.Results.oct_band;
do_plot  = p.Results.plot;
ttl      = p.Results.title;

%% ── Shapes ───────────────────────────────────────────────────────────
J = size(x, 1);        % number of reference sensors
M = size(d, 1);         % number of error mics

fprintf('\n[check_coherence_NR] Computing NR bound for M = %d mic(s), J = %d reference(s)...\n', M, J);

%% ── Pre-compute reference auto/cross spectra Sxx (J x J per freq) ────
% Sxx(:,:,f) is the J x J reference cross-spectral matrix at frequency f.
% Built once; shared by every error mic.
[S_x1x1, f] = pwelch(x(1,:), win, noverlap, nfft, fs);   % to get f and Nf
Nf  = numel(f);
Sxx = zeros(J, J, Nf);
for a = 1:J
    for b = 1:J
        if a == b
            Sxx(a,b,:) = pwelch(x(a,:), win, noverlap, nfft, fs);
        else
            Sxx(a,b,:) = cpsd(x(a,:), x(b,:), win, noverlap, nfft, fs);
        end
    end
end

%% ── Loop over error mics ─────────────────────────────────────────────
for m = 1:M

    d_m = d(m, :);

    % Disturbance auto-spectrum
    S_dd = pwelch(d_m, win, noverlap, nfft, fs);

    if J == 1
        % ---- Single reference: use MATLAB's built-in mscohere directly ----
        % mscohere returns |S_dx|^2/(S_dd*S_xx) = ordinary coherence, exactly
        % the quantity we need. 
        C_dx = mscohere(d_m, x(1,:), win, noverlap, nfft, fs);   % [Nf x 1]
    else
        % ---- Multiple references: multiple coherence (no MATLAB built-in) ----
        %   C(f) = p(f)^H * Sxx(f)^{-1} * p(f) / S_dd(f)
        % Built from cpsd/pwelch spectra. If the J references are highly
        % correlated, Sxx is near-singular and a raw inverse overfits, pinning
        % C to 1 and the NR bound to ~156 dB. Diagonal loading stabilises it.
        p_dx = zeros(J, Nf);
        for j = 1:J
            p_dx(j,:) = cpsd(d_m, x(j,:), win, noverlap, nfft, fs);
        end
        reg  = 1e-6;                         % diagonal-loading factor
        C_dx = zeros(Nf, 1);
        for fi = 1:Nf
            pj   = p_dx(:, fi);                          % [J x 1]
            Sm   = Sxx(:,:,fi);                          % [J x J]
            load = reg * real(trace(Sm)) / J;            % scale to local power
            Sm   = Sm + load * eye(J);                   % diagonal loading
            num  = real(pj' * (Sm \ pj));                % p^H Sxx^-1 p
            C_dx(fi) = num / (S_dd(fi) + eps);
        end
    end
    C_dx = min(max(C_dx, 0), 1);

    % NR bound
    NR_dB = -10*log10(1 - C_dx + eps);

    % Restrict to requested band
    mask    = f >= f_range(1) & f <= f_range(2);
    f_m     = f(mask);
    C_dx_m  = C_dx(mask);
    NR_dB_m = NR_dB(mask);
    S_dd_m  = S_dd(mask);

    results(m).f     = f_m;     %#ok<AGROW>
    results(m).C_dx  = C_dx_m;
    results(m).NR_dB = NR_dB_m;

    % 1/3-octave summary (optional)
    if do_oct
        k_vals = -30:30;
        fc_all = 1000 * 2.^(k_vals/3);
        fc     = fc_all(fc_all >= f_range(1) & fc_all <= f_range(2));
        n_b    = numel(fc);
        C_oct  = NaN(1, n_b);
        NR_oct = NaN(1, n_b);
        for i = 1:n_b
            f_lo = fc(i) * 2^(-1/6);
            f_hi = fc(i) * 2^( 1/6);
            bmask = f_m >= f_lo & f_m <= f_hi;
            if ~any(bmask), continue; end
            w_b = S_dd_m(bmask);  w_b = w_b / (sum(w_b) + eps);
            C_oct(i)  = sum(C_dx_m(bmask) .* w_b);
            NR_oct(i) = -10*log10(1 - C_oct(i) + eps);
        end
        results(m).fc_oct    = fc;
        results(m).C_dx_oct  = C_oct;
        results(m).NR_oct_dB = NR_oct;
    end

    % Summary stats
    results(m).NR_median_dB = median(NR_dB_m(isfinite(NR_dB_m)));
    [best_NR, idx_best]     = max(NR_dB_m);
    results(m).NR_best_dB   = best_NR;
    results(m).NR_best_freq = f_m(idx_best);

    fprintf('  Mic %d: median NR bound = %.1f dB | best = %.1f dB at %.0f Hz\n', ...
            m, results(m).NR_median_dB, best_NR, f_m(idx_best));
end

%% ── Plot ─────────────────────────────────────────────────────────────
if do_plot
    cols = lines(M);
    co_label = sprintf('Coherence  C_{d|x}(\\omega)  (J=%d)', J);

    figure('Name', ['Coherence NR Bound ' ttl], 'Color', 'w');

    subplot(2,1,1); hold on; grid on;
    for m = 1:M
        plot(results(m).f, results(m).C_dx, 'Color', cols(m,:), ...
             'LineWidth', 1.3, 'DisplayName', sprintf('Mic %d', m));
    end
    yline(0.5, 'k:', 'HandleVisibility', 'off');
    ylim([0 1]);
    xlabel('Frequency (Hz)'); ylabel(co_label);
    title(['Reference-Disturbance Coherence   ' ttl]);
    if M > 1, legend('Location','best'); end

    subplot(2,1,2); hold on; grid on;
    for m = 1:M
        plot(results(m).f, results(m).NR_dB, 'Color', cols(m,:), ...
             'LineWidth', 1.3, 'DisplayName', sprintf('Mic %d', m));
    end
    xlabel('Frequency (Hz)'); ylabel('Max achievable NR (dB)');
    title('Coherence-derived NR bound — per-bin');
    if M > 1, legend('Location','best'); end

    if do_oct
        figure('Name', ['Coherence NR Bound — 1-3 Octave ' ttl], 'Color', 'w');
        n_b   = numel(results(1).fc_oct);
        x_pos = 1:n_b;

        subplot(2,1,1); hold on; grid on;
        for m = 1:M
            bar(x_pos + (m-1)*0.8/M - 0.4 + 0.4/M, results(m).C_dx_oct, ...
                0.8/M, 'FaceColor', cols(m,:), 'DisplayName', sprintf('Mic %d', m));
        end
        set(gca, 'XTick', x_pos, 'XTickLabel', round(results(1).fc_oct), ...
                 'XTickLabelRotation', 45, 'FontSize', 8);
        ylim([0 1]);  ylabel('Coherence (1/3-oct avg)');
        title(['1/3-Octave Coherence   ' ttl]);
        if M > 1, legend('Location','best'); end

        subplot(2,1,2); hold on; grid on;
        for m = 1:M
            bar(x_pos + (m-1)*0.8/M - 0.4 + 0.4/M, results(m).NR_oct_dB, ...
                0.8/M, 'FaceColor', cols(m,:), 'DisplayName', sprintf('Mic %d', m));
        end
        set(gca, 'XTick', x_pos, 'XTickLabel', round(results(1).fc_oct), ...
                 'XTickLabelRotation', 45, 'FontSize', 8);
        xlabel('1/3-Octave centre frequency (Hz)');
        ylabel('Max NR bound (dB)');
        title('1/3-Octave NR Bound');
        if M > 1, legend('Location','best'); end
    end
end

end
