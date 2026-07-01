%% =============================================================
%  main.m  —  General J x K x M Active Noise Control
%
%  STRUCTURE:
%    Parameters  ->  Call load/generate  ->  Check coherence/max
%    attenuation ->  Call Mcanc_FxLMS  ->  Call plot results
%
%  This file contains ONLY parameters and function calls. All logic
%  lives in separate function files:
%    load_paths.m, generate_signals.m, check_coherence_NR.m,
%    Mcanc_FxLMS.m, check_anc_status.m, detect_harmonics.m,
%    plot_time_mse.m, plot_psd.m, plot_third_octave.m,
%    plot_spectrogram_anc.m, plot_harmonic_attenuation.m,
%    plot_nr_vs_achieved.m
%
%  DIMENSIONS (fully general J x K x M):
%    J = number of independent REFERENCE sensors   (cfg.J)
%    K = number of control LOUDSPEAKERS            (cfg.K)
%    M = number of error MICROPHONES               (cfg.M)
%
%  TO SWITCH CONFIGURATIONS (1x1x1, 2x2x2, 1x6x6, 2x4x4, etc.):
%    Only Section 1 (System configuration) and Section 2 (path files)
%    need to change. Nothing below Section 2 needs editing.
%
%  uses Sec_hat = S (perfect/clean secondary path knowledge —
%  "Run C" from the wider study) to validate the full pipeline first,
%  before testing imperfect secondary-path identification (I'll do it later).
%% =============================================================

close all; clear; clc;


%% =============================================================
%  SECTION 1 — System configuration  (EDIT THIS TO CHANGE CONFIG)
%% =============================================================

cfg.J  = 3;          % number of reference microphones        (J)
cfg.K  = 4;          % number of control loudspeakers / sec.  (K)
cfg.M  = 2;          % number of error microphones            (M)
cfg.fs = 16000;       % sampling frequency (Hz)
cfg.T  = 60;          % signal duration (s)

cfg.Lw  = 1024;       % FxLMS control filter length (taps)
cfg.muw = 1e-7;       % FxLMS step size (lowered from 1e-6: 3x4x2 with real
                      % measured secondary paths diverged at 1e-6. Tune using
                      % check_paths output: muw < ~2/(Lw*J*M*E_S).)

cfg.f_lo = 25;        % lower frequency bound for analysis / harmonic search
cfg.f_hi = 500;       % upper frequency bound for analysis


%% =============================================================
%  SECTION 2 — Path and signal source configuration
%% =============================================================

cfg.path_source = 'matfile';   % 'matfile' | 'measured' | 'simulated'

% ── BX's data (seminar): 6x6 .mat arrays, error x ref/source x length ──
%   PrimaryPath_6x6.mat   : [6 error x 6 reference x 512]
%   SecondaryPath_6x6.mat : [6 error x 6 source    x 256]
% load_paths auto-detects the stored variable and takes the FIRST
%   M=2 error mics, J=3 references (primary), K=4 sources (secondary).
cfg.primary_mat   = "PrimaryPath_6x6.mat";
cfg.secondary_mat = "SecondaryPath_6x6.mat";

% ── Primary paths P : {M x J} cell — P{m,j} = ref j -> error mic m ─────
% ── Secondary paths S : {M x K} cell — S{m,k} = source k -> mic m ─────
% (only used when cfg.path_source = 'measured')
cfg.path_files.P = {"P_m1_j1.csv","P_m1_j2.csv","P_m1_j3.csv";   % mic 1, refs 1..3
                    "P_m2_j1.csv","P_m2_j2.csv","P_m2_j3.csv"};  % mic 2, refs 1..3

% Secondary paths: {M x K} = {2 x 4} = 8 files (speaker k -> error mic m)
% >>> These are the 8 secondary paths BX measures in the small ANC window <<<
cfg.path_files.S = {"S_m1_k1.csv","S_m1_k2.csv","S_m1_k3.csv","S_m1_k4.csv";   % mic 1, spk 1..4
                    "S_m2_k1.csv","S_m2_k2.csv","S_m2_k3.csv","S_m2_k4.csv"};  % mic 2, spk 1..4

% ── To validate the pipeline WITHOUT measured CSVs, set:
%       cfg.path_source = 'simulated';
%    load_paths.m will then synthesise P:[2x3xLp] and S:[2x4xLs] for you.

% ── Other example configs (uncomment + match Section 1) ───────────────
% 1x1x1 SISO:
% cfg.path_files.P = {"Primary_Path.csv"};
% cfg.path_files.S = {"Spareset left proper calib genset anc test.csv"};

% ── Simulated path settings (used only if cfg.path_source='simulated') ─
cfg.Lp   = 512;
cfg.Ls   = 256;
cfg.p_lo = 80;   cfg.p_hi = 5000;
cfg.s_lo = 80;   cfg.s_hi = 5000;

% ── Reference signal source ─────────────────────────────────────────
cfg.signal_source = 'mixed';   % 'wav' | 'random' | 'synthetic_genset' | 'mixed'

% MIXED: first cfg.n_real references = real wav recording(s); the rest are
% independent synthetic genset-like sources. Here: 1 real Tanglin Halt
% recording + 2 synthetic independent sources = 3 genuinely independent
% references, one of which is the real genset.
cfg.n_real   = 1;       % how many references come from real wav files
cfg.wav_file = { "Tanglin_Halt-Power_Generator_Front(real).wav" };  % the real one(s)

cfg.gen_f0    = [73 41];  % fundamental per SYNTHETIC reference (Hz). With
                          % n_real=1, refs 2 and 3 get 73 Hz and 41 Hz — so the
                          % 3 sources have DISTINCT tones (real genset ~50 Hz,
                          % plus two independent machines at 73 and 41 Hz),
                          % which is physically realistic and visually distinct.
                          % Use a scalar (e.g. 50) to give them all the same f0.
cfg.gen_nharm = 6;      % synthetic harmonics
cfg.gen_bb    = 0.5;    % synthetic broadband floor (fraction of tonal RMS)

% Used only if cfg.signal_source = 'random'
cfg.ref_lo  = 100;
cfg.ref_hi  = 1000;
cfg.ref_snr = 40;

% Optional uncorrelated error-mic noise (dB). With 3 independent sources the
% NR bound is already finite/realistic, so this is OFF by default. Set a
% finite value (e.g. 30) only to model extra error-mic sensor noise.
cfg.dist_snr = Inf;

% Optional: per-reference extra propagation delay (samples) for J>1 'wav'.
% Leave unset to use the default 0,3,6,... Ignored when J=1.
% cfg.ref_delays = [0 4];

rng(1);   % reproducibility


%% =============================================================
%  SECTION 3 — Load paths and generate signals
%% =============================================================

[P, S, Lp, Ls] = load_paths(cfg);          % P:[M x J x Lp], S:[M x K x Ls]

check_paths(P, S, cfg.fs);                  % diagnostics: path gains + muw guide

[x, d] = generate_signals(cfg, P);          % x:[J x N], d:[M x N]


%% =============================================================
%  SECTION 4 — Check coherence / maximum achievable attenuation
%             (BEFORE running any ANC — theoretical ceiling set by
%              sensor placement; multiple coherence when J>1, per
%              Kuo & Morgan Eq. 3.1.4 / 3.1.11)
%% =============================================================

coherence_results = check_coherence_NR(x, d, cfg.fs, ...
                                        'f_range', [cfg.f_lo cfg.f_hi], ...
                                        'title', sprintf('J=%d K=%d M=%d', cfg.J, cfg.K, cfg.M));


%% =============================================================
%  SECTION 5 — Run Mcanc_FxLMS  (Sec_hat = S : perfect knowledge,
%             validating the pipeline before testing imperfect ID)
%% =============================================================

Sec_hat = S;   % perfect / clean secondary path estimate ("Run C")

[e, Wc, n_div, y, yp] = Mcanc_FxLMS(cfg, x, d, S, Sec_hat);


%% =============================================================
%  SECTION 5b — Dimension check  (shows the correct variable counts;
%               this is what the seminar demo is meant to verify)
%% =============================================================

fprintf('\n=== Variable dimension check (J=%d, K=%d, M=%d, Lw=%d) ===\n', ...
        cfg.J, cfg.K, cfg.M, cfg.Lw);
fprintf('  x  (reference)        : [%d x %d]        (J x N)\n',     size(x,1),  size(x,2));
fprintf('  d  (disturbance)      : [%d x %d]        (M x N)\n',     size(d,1),  size(d,2));
fprintf('  P  (primary paths)    : [%d x %d x %d]   (M x J x Lp) -> %d primary paths\n', ...
        size(P,1), size(P,2), size(P,3), cfg.M*cfg.J);
fprintf('  S  (secondary paths)  : [%d x %d x %d]   (M x K x Ls) -> %d secondary paths\n', ...
        size(S,1), size(S,2), size(S,3), cfg.M*cfg.K);
fprintf('  Wc (control filters)  : [%d x %d x %d]   (K x J x Lw) -> %d control filters, %d taps each\n', ...
        size(Wc,1), size(Wc,2), size(Wc,3), cfg.K*cfg.J, cfg.Lw);
fprintf('  y  (speaker drive)    : [%d x %d]        (K x N)\n',     size(y,1),  size(y,2));
fprintf('  e  (error)            : [%d x %d]        (M x N)\n',     size(e,1),  size(e,2));
fprintf('  yp (anti-noise)       : [%d x %d]        (M x N)\n',     size(yp,1), size(yp,2));
fprintf('===========================================================\n');


%% =============================================================
%  SECTION 5c — Automatic variable-count verification
%   Fails loudly (assert) if ANY array does not match the J,K,M set
%   in Section 1. Pure check — no effect on results. If you see the
%   "[verify] All variable counts CORRECT" line, everything matched.
%% =============================================================

assert(isequal(size(x),  [cfg.J cfg.fs*cfg.T]),  'x must be [J x N] = [%d x %d]', cfg.J, cfg.fs*cfg.T);
assert(isequal(size(d),  [cfg.M cfg.fs*cfg.T]),  'd must be [M x N] = [%d x %d]', cfg.M, cfg.fs*cfg.T);
assert(size(P,1)==cfg.M && size(P,2)==cfg.J,     'P must be [M x J x Lp] = [%d x %d x .]', cfg.M, cfg.J);
assert(size(S,1)==cfg.M && size(S,2)==cfg.K,     'S must be [M x K x Ls] = [%d x %d x .]', cfg.M, cfg.K);
assert(isequal(size(Wc),[cfg.K cfg.J cfg.Lw]),   'Wc must be [K x J x Lw] = [%d x %d x %d]', cfg.K, cfg.J, cfg.Lw);
assert(isequal(size(y),  [cfg.K cfg.fs*cfg.T]),  'y must be [K x N] = [%d x %d]', cfg.K, cfg.fs*cfg.T);
assert(isequal(size(e),  [cfg.M cfg.fs*cfg.T]),  'e must be [M x N] = [%d x %d]', cfg.M, cfg.fs*cfg.T);
assert(isequal(size(yp), [cfg.M cfg.fs*cfg.T]),  'yp must be [M x N] = [%d x %d]', cfg.M, cfg.fs*cfg.T);

fprintf('\n[verify] All variable counts CORRECT for J=%d, K=%d, M=%d:\n', cfg.J, cfg.K, cfg.M);
fprintf('   control filters = K*J = %d  (each %d taps)\n', cfg.K*cfg.J, cfg.Lw);
fprintf('   secondary paths = M*K = %d\n', cfg.M*cfg.K);
fprintf('   primary paths   = M*J = %d\n', cfg.M*cfg.J);
fprintf('[verify] Every array matched its J/K/M rule.\n');


%% =============================================================
%  SECTION 6 — Plot results  (time, PSD, spectrum, 1/3-oct,
%             spectrogram, harmonic bars, coherence-vs-achieved)
%% =============================================================

run_label = sprintf('J=%d K=%d M=%d (perfect S_hat)', cfg.J, cfg.K, cfg.M);

idx_ss            = check_anc_status(cfg, d, e, n_div, 'label', run_label);
[f0, harm_freqs]  = detect_harmonics(cfg, x);

plot_time_mse(cfg, e, n_div, 'label', run_label);

plot_psd(cfg, d, e, idx_ss, 'label', run_label, 'harm_freqs', harm_freqs);

for m = 1:cfg.M
    plot_third_octave(d, e, cfg.fs, 'mic', m, 'f_range', [cfg.f_lo cfg.f_hi], ...
                       'title', sprintf('%s Mic %d', run_label, m));
end

for m = 1:cfg.M
    plot_spectrogram_anc(d, e, cfg.fs, 'mic', m, 'f_max', cfg.f_hi, ...
                         'title', sprintf('%s Mic %d', run_label, m));
end

plot_harmonic_attenuation(d, e, cfg.fs, f0, 'mic', 1, 'n_harm', 6, ...
                           'title', [run_label ' Mic 1']);

plot_nr_vs_achieved(cfg, d, e, idx_ss, coherence_results, ...
                     'label', run_label, 'harm_freqs', harm_freqs);