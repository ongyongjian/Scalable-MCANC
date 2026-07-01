function idx_ss = check_anc_status(cfg, d, e, n_div, varargin)
%% check_anc_status  Steady-state window + stable/diverged attenuation summary
%
%  WHAT IT DOES
%  ------------
%  Determines the steady-state sample window (last 5 s if stable, or up
%  to divergence if not), prints whether the run was stable or diverged,
%  and if stable, prints the achieved attenuation per error mic.
%
%  INPUTS
%  ------
%  cfg      : struct — needs fs, T, M
%  d        : [M x N]  disturbance
%  e        : [M x N]  error signal
%  n_div    : scalar — divergence sample index, NaN if stable
%  varargin : optional name-value pairs
%             'label' , char — run label for printed messages
%
%  OUTPUT
%  ------
%  idx_ss : steady-state sample index range, for use by downstream plots
%
%  CALL EXAMPLE
%  ------------
%  idx_ss = check_anc_status(cfg, d, e, n_div, 'label', 'Run C');

p = inputParser;
addParameter(p, 'label', 'ANC Run', @ischar);
parse(p, varargin{:});
lbl = p.Results.label;

N = cfg.fs * cfg.T;

if isnan(n_div)
    idx_ss = (N - 5*cfg.fs + 1) : N;
    fprintf('\n[check_anc_status] %s: stable run.\n', lbl);

    for m = 1:cfg.M
        atten_m = 20*log10(rms(d(m,idx_ss)) / (rms(e(m,idx_ss)) + eps));
        fprintf('  Mic %d attenuation: %+.2f dB\n', m, atten_m);
    end
else
    idx_ss = 1 : max(2, n_div-1);
    fprintf('\n[check_anc_status] %s: DIVERGED at t=%.2fs — using pre-divergence window only.\n', ...
            lbl, n_div/cfg.fs);
end

end