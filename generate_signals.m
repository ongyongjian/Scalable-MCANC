function [x, d] = generate_signals(cfg, P)
%% generate_signals  Build reference x [J x N] and disturbance d [M x N].
%
%  General J x K x M:
%    J = number of independent reference sensors (rows of x)
%    M = number of error microphones (rows of d)
%  Reduces to the old single-reference behaviour when J = 1.
%
%  The disturbance at mic m is the sum over all references of each
%  reference filtered through its primary path:
%     d(m,:) = Σ_j filter(P(m,j,:), 1, x(j,:))
%  Does not depend on K (loudspeaker count) at all.
%
%  INPUTS
%  ------
%  cfg.signal_source : 'wav' or 'random'
%  cfg.wav_file       : for 'wav' — either a single filename (J=1) or a
%                       {1 x J} cell of filenames (one independent source
%                       recording per reference) when J>1
%  cfg.ref_lo/hi/snr   : used only if signal_source='random'
%  cfg.fs, cfg.T        : sampling rate and duration
%  cfg.J                : number of reference sensors (optional; default 1)
%  cfg.dist_snr         : optional error-mic noise (dB); default Inf (off)
%  P                    : [M x J x Lp] primary path impulse responses
%
%  OUTPUTS
%  -------
%  x : [J x N]  reference signals (one row per reference source)
%  d : [M x N]  disturbance at each error mic

    N = cfg.fs * cfg.T;

    if isfield(cfg, 'J'), J = cfg.J; else, J = 1; end

    % P may arrive as [M x J x Lp] (general) or [M x Lp] (legacy J=1).
    if ndims(P) == 2
        P = reshape(P, size(P,1), 1, size(P,2));   % [M x 1 x Lp]
    end
    assert(size(P,2) == J, ...
        'generate_signals: P second dim (%d) must equal J (%d).', size(P,2), J);

    x = zeros(J, N);

    switch cfg.signal_source
        case 'wav'
            % J independent reference sources, one WAV file per reference.
            % cfg.wav_file may be:
            %   - a string/char (J=1): single file, OR
            %   - a {1 x J} or {J x 1} cell of filenames (J>1): one per ref.
            % Each file is an independent source recording; reference j is
            % that recording, and it reaches the error mics through its own
            % primary path P(:,j,:). Because the J sources are genuinely
            % independent, the J references carry J real degrees of freedom
            % (no singular Sxx, no 156 dB artefact).
            wf = cfg.wav_file;
            if ischar(wf) || isstring(wf), wf = {wf}; end   % wrap single file
            assert(numel(wf) == J, ...
                ['generate_signals: need J=%d wav file(s) in cfg.wav_file, ' ...
                 'got %d. Provide one recording per reference source.'], J, numel(wf));

            for j = 1:J
                [x_raw, fs_in] = audioread(wf{j});
                if size(x_raw, 2) > 1, x_raw = mean(x_raw, 2); end   % to mono
                x_res = resample(x_raw, cfg.fs, fs_in);
                if length(x_res) < N
                    error('generate_signals: WAV "%s" shorter than cfg.T = %d s.', ...
                          wf{j}, cfg.T);
                end
                x(j,:) = x_res(1:N).';
            end

        case 'random'
            % Independent band-limited reference per sensor.
            b = fir1(63, [cfg.ref_lo cfg.ref_hi] / (cfg.fs/2));
            for j = 1:J
                x_clean = filter(b, 1, randn(N, 1));
                noise   = randn(N, 1);
                noise   = noise / (rms(noise) + eps) * rms(x_clean) / 10^(cfg.ref_snr/20);
                x(j,:)  = (x_clean + noise).';
            end

        case 'mixed'
            % First cfg.n_real references come from real wav recording(s);
            % the remaining (J - n_real) are independent synthetic genset-like
            % sources. Use when you have one real genset recording but want a
            % genuine multi-reference demo. Real and synthetic references are
            % independent, and the synthetic broadband keeps the set stable.
            if isfield(cfg,'n_real'), n_real = cfg.n_real; else, n_real = 1; end
            wf = cfg.wav_file; if ischar(wf)||isstring(wf), wf = {wf}; end
            assert(numel(wf) >= n_real, ...
                'generate_signals: need >= n_real=%d wav file(s) in cfg.wav_file.', n_real);
            % --- real references from wav ---
            for j = 1:n_real
                [w, fin] = audioread(wf{j});
                if size(w,2) > 1, w = mean(w,2); end
                w = resample(w, cfg.fs, fin);
                assert(numel(w) >= N, 'generate_signals: wav "%s" shorter than cfg.T.', wf{j});
                x(j,:) = w(1:N).';
            end
            % --- remaining references: independent synthetic genset-like ---
            if isfield(cfg,'gen_f0'),    f0v=cfg.gen_f0;      else, f0v=50;    end
            if isfield(cfg,'gen_nharm'), nharm=cfg.gen_nharm; else, nharm=6;   end
            if isfield(cfg,'gen_bb'),    bb_frac=cfg.gen_bb;  else, bb_frac=0.5; end
            tt = (0:N-1)/cfg.fs;  bbb = fir1(63,[40 2000]/(cfg.fs/2));
            for j = n_real+1 : J
                % per-reference fundamental: f0v may be scalar (shared) or a
                % vector with one fundamental per SYNTHETIC reference, so each
                % independent source can have its own tones (e.g. 50, 73, 41).
                if isscalar(f0v), f0 = f0v; else, f0 = f0v(j - n_real); end
                sig = zeros(1,N);
                for h = 1:nharm
                    amp = (1/h)*(0.7+0.6*rand);
                    ph  = 2*pi*rand + 0.5*cumsum(randn(1,N))/cfg.fs;
                    sig = sig + amp*sin(2*pi*f0*h*tt + ph);
                end
                bb = filter(bbb,1,randn(1,N));  bb = bb/(rms(bb)+eps)*rms(sig)*bb_frac;
                x(j,:) = sig + bb;
                x(j,:) = x(j,:)/(rms(x(j,:))+eps)*rms(x(1,:));   % match real level
            end

        case 'synthetic_genset'
            % J INDEPENDENT genset-like reference sources. Each reference has
            % the SAME tonal STRUCTURE (fundamental cfg.gen_f0 + harmonics) but
            % its OWN independent amplitudes, slow random phase drift, and its
            % OWN broadband floor. Result: references are genset-like in
            % spectrum yet genuinely independent (low mutual coherence,
            % well-conditioned Sxx), so the J>1 architecture is functionally
            % meaningful. The broadband floor also spreads the power, which
            % keeps FxLMS stable (a pure tone concentrates power and diverges).
            %
            % Use this when only single-source field recordings are available
            % but you want to demonstrate genuine multi-reference operation.
            if isfield(cfg,'gen_f0'),    f0v   = cfg.gen_f0;    else, f0v   = 50;  end
            if isfield(cfg,'gen_nharm'), nharm = cfg.gen_nharm; else, nharm = 6;   end
            if isfield(cfg,'gen_bb'),    bb_frac = cfg.gen_bb;  else, bb_frac = 0.5; end
            tt   = (0:N-1) / cfg.fs;
            bbb  = fir1(63, [40 2000] / (cfg.fs/2));
            for j = 1:J
                % per-reference fundamental (scalar = shared, or vector with
                % one fundamental per reference so each source has its own tones)
                if isscalar(f0v), f0 = f0v; else, f0 = f0v(j); end
                sig = zeros(1, N);
                for h = 1:nharm
                    amp = (1/h) * (0.7 + 0.6*rand);                  % independent amplitude
                    ph  = 2*pi*rand + 0.5*cumsum(randn(1,N))/cfg.fs; % independent slow phase drift
                    sig = sig + amp*sin(2*pi*f0*h*tt + ph);
                end
                bb  = filter(bbb, 1, randn(1, N));                   % independent broadband
                bb  = bb / (rms(bb)+eps) * rms(sig) * bb_frac;
                x(j,:) = sig + bb;
            end
            x = x / max(abs(x(:))) * 0.5;   % normalise headroom

        otherwise
            error('generate_signals: unknown cfg.signal_source ''%s''', cfg.signal_source);
    end

    %% ── Disturbance: sum each reference through its primary path ──────
    d = zeros(cfg.M, N);
    for m = 1:cfg.M
        for j = 1:J
            p_mj   = reshape(P(m,j,:), 1, []);     % [1 x Lp]
            d(m,:) = d(m,:) + filter(p_mj, 1, x(j,:));
        end
    end

    %% ── Optional uncorrelated error-mic noise ────────────────────────
    % With J independent sources the references already give a finite,
    % realistic coherence bound, so this is OFF by default (Inf). Set
    % cfg.dist_snr to a finite dB value only if you want to model extra
    % sensor noise at the error mics on top of that.
    if isfield(cfg, 'dist_snr'), dist_snr = cfg.dist_snr; else, dist_snr = Inf; end
    if isfinite(dist_snr)
        for m = 1:cfg.M
            w      = randn(1, N);
            w      = w / (rms(w) + eps) * rms(d(m,:)) / 10^(dist_snr/20);
            d(m,:) = d(m,:) + w;
        end
    end

    fprintf('\n--- Signal dimensions ---\n');
    fprintf('  x : [%d x %d]   (J x N)\n', size(x,1), size(x,2));
    fprintf('  d : [%d x %d]   (M x N)\n', size(d,1), size(d,2));
    if isfinite(dist_snr)
        fprintf('  (error-mic disturbance SNR = %.0f dB; coherence bound will be finite)\n', dist_snr);
    else
        fprintf('  (dist_snr = Inf: deterministic disturbance, coherence = 1)\n');
    end

    %% ---------- reference signal plots (time + frequency) ----------
    Jr = size(x, 1);
    leg = strings(1, Jr);
    for j = 1:Jr, leg(j) = sprintf('ref %d', j); end

    % ===== reference signals — time domain (first 0.1 s, overlaid) =====
    figure('Name','Reference signals — time domain','Color','w'); hold on;
    nshow = min(N, round(0.1 * cfg.fs));        % first 100 ms (full N is too dense)
    tt = (0:nshow-1) / cfg.fs * 1000;            % ms
    for j = 1:Jr
        plot(tt, x(j, 1:nshow), 'LineWidth', 1);
    end
    hold off; grid on; xlabel('Time (ms)'); ylabel('Amplitude');
    title('Reference signals — time domain (first 100 ms)');
    legend(leg, 'Location', 'best'); legend boxoff;

    % ===== reference signals — frequency domain (PSD, overlaid) =====
    figure('Name','Reference signals — frequency domain','Color','w'); hold on;
    nfft = 8192; win = hann(4096); nov = 2048;
    for j = 1:Jr
        [Pxx, fpx] = pwelch(x(j,:), win, nov, nfft, cfg.fs);
        plot(fpx, 10*log10(Pxx + eps), 'LineWidth', 1);
    end
    hold off; grid on; xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
    title('Reference signals — frequency domain');
    legend(leg, 'Location', 'best'); legend boxoff;
    xlim([0 min(cfg.fs/2, 2000)]);   % show tones + broadband up to 2 kHz
end