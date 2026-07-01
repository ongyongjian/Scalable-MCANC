function check_paths(P, S, fs)
%% check_paths  Inspect primary (P) and secondary (S) paths: prints gains,
%               delays and the FxLMS step-size guide, and plots the time and
%               frequency responses (all paths overlaid on one axes each).
%
%  USAGE (after load_paths in main.m):
%     check_paths(P, S, cfg.fs)
%
%  P : [M x J x Lp] primary paths   (P(m,j,:) = ref j -> error mic m)
%  S : [M x K x Ls] secondary paths (S(m,k,:) = source k -> error mic m)
%
%  FIGURES PRODUCED
%    1) all primary   impulse responses (time)      overlaid
%    2) all primary   frequency responses (mag, dB) overlaid
%    3) all secondary impulse responses (time)      overlaid
%    4) all secondary frequency responses (mag, dB) overlaid

    [M, J, Lp] = size(P);
    [~, K, Ls] = size(S);

    %% ---------- printed diagnostics ----------
    fprintf('\n================ PATH DIAGNOSTICS ================\n');
    fprintf('P: [M=%d x J=%d x Lp=%d]   S: [M=%d x K=%d x Ls=%d]\n', M,J,Lp, M,K,Ls);
    if any(~isfinite(P(:))), fprintf('*** WARNING: P contains NaN/Inf! ***\n'); end
    if any(~isfinite(S(:))), fprintf('*** WARNING: S contains NaN/Inf! ***\n'); end

    fprintf('\n--- PRIMARY paths P(m,j) ---\n');
    fprintf('   %-8s %-10s %-10s %-10s %-8s\n','(m,j)','peak','rms','energy','delay');
    for m = 1:M
        for j = 1:J
            h = squeeze(P(m,j,:));
            [pk, ipk] = max(abs(h));
            fprintf('   (%d,%d)    %-10.4g %-10.4g %-10.4g %-8d\n', m,j, pk, rms(h), sum(h.^2), ipk-1);
        end
    end

    fprintf('\n--- SECONDARY paths S(m,k)  [set the FxLMS stability limit] ---\n');
    fprintf('   %-8s %-10s %-10s %-10s %-8s\n','(m,k)','peak','rms','energy','delay');
    S_energy_max = 0;
    for m = 1:M
        for k = 1:K
            h = squeeze(S(m,k,:));
            [pk, ipk] = max(abs(h));
            e = sum(h.^2);  S_energy_max = max(S_energy_max, e);
            fprintf('   (%d,%d)    %-10.4g %-10.4g %-10.4g %-8d\n', m,k, pk, rms(h), e, ipk-1);
        end
    end

    fprintf('\n--- Stability guide ---\n');
    fprintf('   Largest secondary-path energy = %.4g\n', S_energy_max);
    Lw = 1024;
    mu_bound = 2 / (Lw * J * M * S_energy_max + eps);
    fprintf('   For Lw=%d, J=%d, M=%d:  muw should be roughly < %.2e\n', Lw,J,M, mu_bound);
    fprintf('   (use ~10x smaller for safety:  muw ~ %.2e)\n', mu_bound/10);
    fprintf('==================================================\n\n');

    %% ---------- plots ----------
    nfft = max(2048, 2^nextpow2(max(Lp,Ls)));
    f    = (0:nfft/2-1) * fs/nfft;
    tP   = (0:Lp-1) / fs * 1000;    % ms
    tS   = (0:Ls-1) / fs * 1000;    % ms

    % ---- 1) PRIMARY impulse responses (time) ----
    figure('Name','Primary paths — time'); hold on;
    legP = cell(1, M*J); c = 0;
    for m = 1:M
        for j = 1:J
            c = c + 1;
            plot(tP, squeeze(P(m,j,:)), 'LineWidth', 1);
            legP{c} = sprintf('P(m%d,j%d)', m, j);
        end
    end
    grid on; xlabel('Time (ms)'); ylabel('Amplitude');
    title('Primary path impulse responses P(m,j)'); legend(legP,'Location','best');

    % ---- 2) PRIMARY frequency responses (magnitude dB) ----
    figure('Name','Primary paths — frequency'); hold on;
    for m = 1:M
        for j = 1:J
            H = fft(squeeze(P(m,j,:)), nfft);
            plot(f, 20*log10(abs(H(1:nfft/2)) + eps), 'LineWidth', 1);
        end
    end
    grid on; xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    title('Primary path frequency responses P(m,j)'); legend(legP,'Location','best');
    xlim([0 fs/2]);

    % ---- 3) SECONDARY impulse responses (time) ----
    figure('Name','Secondary paths — time'); hold on;
    legS = cell(1, M*K); c = 0;
    for m = 1:M
        for k = 1:K
            c = c + 1;
            plot(tS, squeeze(S(m,k,:)), 'LineWidth', 1);
            legS{c} = sprintf('S(m%d,k%d)', m, k);
        end
    end
    grid on; xlabel('Time (ms)'); ylabel('Amplitude');
    title('Secondary path impulse responses S(m,k)'); legend(legS,'Location','best');

    % ---- 4) SECONDARY frequency responses (magnitude dB) ----
    figure('Name','Secondary paths — frequency'); hold on;
    for m = 1:M
        for k = 1:K
            H = fft(squeeze(S(m,k,:)), nfft);
            plot(f, 20*log10(abs(H(1:nfft/2)) + eps), 'LineWidth', 1);
        end
    end
    grid on; xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    title('Secondary path frequency responses S(m,k)'); legend(legS,'Location','best');
    xlim([0 fs/2]);
end