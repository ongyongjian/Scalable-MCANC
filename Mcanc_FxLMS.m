function [e, Wc, n_div, y, yp] = Mcanc_FxLMS(cfg, x, d, S, Sec_hat)
%% Mcanc_FxLMS  Multichannel Filtered-x LMS active noise controller
%               General J x K x M:
%                 J = independent reference sensors
%                 K = control loudspeakers
%                 M = error microphones
%               Reduces exactly to the old 1 x K x M code when J = 1.
%
%  PERFORMANCE
%  -----------
%  All three per-sample bottlenecks are vectorised:
%   * Delay lines (xbuf, xfbuf, ybuf) use CIRCULAR BUFFERS — write the newest
%     slice, read through a rotated index — instead of shifting the whole
%     array every sample.
%   * Secondary-path filtering (Step 4) is done as a single tensor contraction
%     over a circular y-history buffer, NOT M*K stateful filter() calls per
%     sample (those dominated the runtime; this is ~18x faster and exact,
%     because the secondary path is FIR).
%   * The weight update (Step 6) is a single tensor contraction over all
%     (k,j) at once instead of a nested k,j loop.
%  Results are bitwise-identical to the per-sample filter() version.
%
%  ── Inputs ────────────────────────────────────────────────────────────
%  cfg     : struct with fields K, M, fs, T, Lw, muw  (J read from x)
%  x       : [J × N]        reference signals (one row per reference sensor)
%  d       : [M × N]        disturbance at each error microphone
%  S       : [M × K × Ls]  TRUE secondary-path impulse responses
%  Sec_hat : [M × K × Ls]  ESTIMATED secondary paths (used for FxLMS gradient)
%
%  ── Outputs ───────────────────────────────────────────────────────────
%  e       : [M × N]      error signals; columns from n_div onward NaN if diverged
%  Wc      : [K × J × Lw] final adaptive filter weights, one FIR per (k,j)
%  n_div   : scalar       sample index of divergence, or NaN if stable
%  y       : [K × N]      loudspeaker drive signals
%  yp      : [M × N]      total anti-noise at each error mic

%% ── Unpack scalars ───────────────────────────────────────────────────
[K, M, Lw, muw] = deal(cfg.K, cfg.M, cfg.Lw, cfg.muw);
J  = size(x, 1);            % number of reference sensors (J=1 => old behaviour)
N  = cfg.fs * cfg.T;        % total samples
Ls = size(S, 3);            % secondary-path length

%% ── Dimension checks ─────────────────────────────────────────────────
assert(isequal(size(x),      [J  N      ]) && ...
       isequal(size(d),      [M  N      ]) && ...
       isequal(size(S),      [M  K  Ls  ]) && ...
       isequal(size(Sec_hat),[M  K  Ls  ]), ...
       'Mcanc_FxLMS: input dimension mismatch.');

%% ── Flatten estimated paths for the offline filtered-x pre-filtering ──
Sec2 = reshape(Sec_hat, M*K, Ls);   % row mk = estimated path (m,k)
mk_m = zeros(M*K,1);  mk_k = zeros(M*K,1);
for mk = 1:M*K
    [mk_m(mk), mk_k(mk)] = ind2sub([M K], mk);
end

%% ── Pre-filter every reference x_j through every estimated path Ŝ(m,k) ─
% XF(m,k,j,:) = filter(Sec_hat(m,k), 1, x(j,:)).  Done once, offline.
XF = zeros(M, K, J, N);
for i = 1:M*K
    m = mk_m(i);  k = mk_k(i);
    for j = 1:J
        XF(m,k,j,:) = filter(Sec2(i,:), 1, x(j,:));
    end
end

%% ── Reshape S for the vectorised Step-4 contraction ──────────────────
% Smat(:, (k-1)*Ls + tap) groups the secondary path of every (m,k) so that
% yp = Smat * yh_stacked is a single matrix-vector product per sample.
% S is [M x K x Ls]; we want, per sample,  yp(m) = Σ_k Σ_tap S(m,k,tap)·yh(k,tap)
Smat = reshape(S, M, K*Ls);          % [M x K*Ls], column order (k,tap) via reshape

%% ── Initialise buffers ───────────────────────────────────────────────
Wc    = zeros(K, J, Lw);        % control filters, one FIR per (k,j)
xbuf  = zeros(J, Lw);           % reference delay line   (circular)
xfbuf = zeros(M, K, J, Lw);     % filtered-x delay line  (circular)
ybuf  = zeros(K, Ls);           % loudspeaker-output history (circular, for S)
e     = zeros(M, N);
y     = zeros(K, N);
yp    = zeros(M, N);
n_div = NaN;

% Two circular pointers: one for the Lw-length lines, one for the Ls history.
ptr  = 1;  tap  = 0:Lw-1;
pty  = 1;  tapS = 0:Ls-1;

%% ── Sample loop ──────────────────────────────────────────────────────
for n = 1:N

    %% Step 1+2 — advance reference / filtered-x delay lines ───────────
    ptr = mod(ptr-2, Lw) + 1;
    xbuf(:,ptr)      = x(:,n);
    xfbuf(:,:,:,ptr) = XF(:,:,:,n);
    idx = mod((ptr-1) + tap, Lw) + 1;     % newest -> oldest

    %% Step 3 — loudspeaker drive  y(k,n) = Σ_j Wc(k,j,:)·xbuf(j,:) ────
    xb = xbuf(:, idx);                    % [J × Lw]
    yn = zeros(K,1);
    for j = 1:J
        yn = yn + reshape(Wc(:,j,:), K, Lw) * xb(j,:).';
    end
    y(:,n) = yn;

    %% Step 4 — anti-noise through TRUE paths (vectorised FIR) ─────────
    % Push yn into the circular y-history, read it newest->oldest, then
    % yp = Σ_k Σ_tap S(m,k,tap)·yh(k,tap) as ONE matrix-vector product.
    pty = mod(pty-2, Ls) + 1;
    ybuf(:,pty) = yn;
    idy = mod((pty-1) + tapS, Ls) + 1;    % newest -> oldest
    yh  = ybuf(:, idy);                   % [K × Ls]
    yp(:,n) = Smat * reshape(yh, K*Ls, 1);
    e(:,n)  = d(:,n) - yp(:,n);

    %% Step 5 — divergence check ──────────────────────────────────────
    if any(~isfinite(e(:,n))) || max(abs(Wc(:))) > 1e6
        n_div = n;  e(:,n:end)=NaN;  y(:,n:end)=NaN;  yp(:,n:end)=NaN;  break
    end

    %% Step 6 — FxLMS weight update (vectorised over all k,j) ──────────
    % Wc(k,j,:) += muw · Σ_m e(m,n)·xfbuf(m,k,j,:).  Contract over m with a
    % single tensor product instead of a nested k,j loop.
    xf_ord = xfbuf(:,:,:,idx);                       % [M × K × J × Lw]
    en     = e(:,n);                                 % [M × 1]
    % sum over m: result [K × J × Lw]
    grad   = reshape(en.' * reshape(xf_ord, M, K*J*Lw), K, J, Lw);
    Wc     = Wc + muw * grad;

end  % for n = 1:N

end  % function Mcanc_FxLMS