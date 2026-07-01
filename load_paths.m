function [P, S, Lp, Ls] = load_paths(cfg)
%% load_paths  Load or simulate primary (P) and secondary (S) paths.
%
%  General J x K x M:
%    J = number of independent reference sensors
%    K = number of control loudspeakers (secondary sources)
%    M = number of error microphones
%
%  THREE SOURCES (cfg.path_source):
%    'matfile'   : load BX-style .mat files holding 3-D arrays
%                  (error x reference/source x length) and slice the first
%                  M / K / J channels. THIS is the seminar path.
%    'measured'  : read individual CSV files via {M x J} / {M x K} cells.
%    'simulated' : synthesise FIR paths (no data files needed).
%
%  OUTPUTS
%  -------
%  P  : [M x J x Lp]  primary path impulse responses (ref j -> mic m)
%  S  : [M x K x Ls]  secondary path impulse responses (source k -> mic m)
%  Lp, Ls : scalar filter lengths read from the loaded data

    if isfield(cfg, 'J'), J = cfg.J; else, J = 1; end
    M = cfg.M;  K = cfg.K;

    switch cfg.path_source

        % =========================================================
        case 'matfile'
        % =========================================================
        % BX's data:
        %   PrimaryPath_6x6.mat   : [6 error x 6 reference x 512]
        %   SecondaryPath_6x6.mat : [6 error x 6 source    x 256]
        % We auto-detect the stored variable, confirm the axis order,
        % then take the first M error mics, first J references (primary)
        % and first K sources (secondary).

            Praw = local_load_3d(cfg.primary_mat,   'primary');
            Sraw = local_load_3d(cfg.secondary_mat, 'secondary');

            % --- report what was found, so axis order can be eyeballed ---
            fprintf('\n--- Loaded .mat arrays (raw, before slicing) ---\n');
            fprintf('  primary   "%s" : [%d x %d x %d]  (expect error x reference x length)\n', ...
                    cfg.primary_mat,   size(Praw,1), size(Praw,2), size(Praw,3));
            fprintf('  secondary "%s" : [%d x %d x %d]  (expect error x source x length)\n', ...
                    cfg.secondary_mat, size(Sraw,1), size(Sraw,2), size(Sraw,3));

            % --- sanity checks on available channels ---
            assert(size(Praw,1) >= M && size(Praw,2) >= J, ...
                ['load_paths: primary array is [%d x %d x .]; need at least ' ...
                 'M=%d error x J=%d reference. Check axis order / file.'], ...
                size(Praw,1), size(Praw,2), M, J);
            assert(size(Sraw,1) >= M && size(Sraw,2) >= K, ...
                ['load_paths: secondary array is [%d x %d x .]; need at least ' ...
                 'M=%d error x K=%d source. Check axis order / file.'], ...
                size(Sraw,1), size(Sraw,2), M, K);

            % --- slice first M / J / K channels ---
            P  = Praw(1:M, 1:J, :);     % [M x J x Lp]
            S  = Sraw(1:M, 1:K, :);     % [M x K x Ls]
            Lp = size(P, 3);
            Ls = size(S, 3);

        % =========================================================
        case 'measured'
        % =========================================================
            PF = cfg.path_files.P;
            if J == 1 && (isrow(PF) || iscolumn(PF)) && numel(PF) == M
                PF = reshape(PF, M, 1);   % legacy {1 x M} -> {M x 1}
            end
            assert(isequal(size(PF), [M, J]), ...
                'load_paths: path_files.P must be a {M x J} = {%d x %d} cell.', M, J);

            p_first = readmatrix(PF{1,1});  p_first = p_first(:).';
            Lp      = length(p_first);
            P       = zeros(M, J, Lp);
            P(1,1,:) = p_first;
            for m = 1:M
                for j = 1:J
                    if m == 1 && j == 1, continue; end
                    p_mj = readmatrix(PF{m,j});  p_mj = p_mj(:).';
                    if length(p_mj) ~= Lp
                        error('load_paths: P{%d,%d} length (%d) ~= P{1,1} (%d).', ...
                            m, j, length(p_mj), Lp);
                    end
                    P(m,j,:) = p_mj;
                end
            end

            s_first  = readmatrix(cfg.path_files.S{1,1});  s_first = s_first(:);
            Ls       = length(s_first);
            S        = zeros(M, K, Ls);
            S(1,1,:) = s_first;
            for m = 1:M
                for k = 1:K
                    if m == 1 && k == 1, continue; end
                    s_mk = readmatrix(cfg.path_files.S{m,k});  s_mk = s_mk(:);
                    if length(s_mk) ~= Ls
                        error('load_paths: S{%d,%d} length (%d) ~= S{1,1} (%d).', ...
                            m, k, length(s_mk), Ls);
                    end
                    S(m,k,:) = s_mk;
                end
            end

        % =========================================================
        case 'simulated'
        % =========================================================
            Lp = cfg.Lp;   Ls = cfg.Ls;
            P  = zeros(M, J, Lp);
            S  = zeros(M, K, Ls);
            bP = fir1(Lp-1, [cfg.p_lo cfg.p_hi] / (cfg.fs/2));
            bS = fir1(Ls-1, [cfg.s_lo cfg.s_hi] / (cfg.fs/2));
            for m = 1:M
                for j = 1:J
                    delay_p = 5*(m-1) + 3*(j-1);
                    gain_p  = (1 - 0.08*(m-1)) * (1 - 0.05*(j-1));
                    P(m,j,:) = gain_p * [zeros(1, delay_p), bP(1:Lp - delay_p)];
                end
            end
            for m = 1:M
                for k = 1:K
                    delay_s  = 7*(m-1) + 11*(k-1);
                    gain_s   = 1 / (1 + 0.20*abs(m - k));
                    s_mk     = gain_s * [zeros(1, delay_s), bS(1:Ls - delay_s)];
                    s_mk     = s_mk + 0.005 * randn(size(s_mk));
                    S(m,k,:) = s_mk;
                end
            end

        otherwise
            error('load_paths: unknown cfg.path_source ''%s''', cfg.path_source);
    end

    %% ── Report final (sliced) dimensions ─────────────────────────────
    fprintf('\n--- Path dimensions (after selection) ---\n');
    fprintf('  P : [%d x %d x %d]  (M x J x Lp)  -> %d primary paths\n', ...
            size(P,1), size(P,2), size(P,3), M*J);
    fprintf('  S : [%d x %d x %d]  (M x K x Ls)  -> %d secondary paths\n', ...
            size(S,1), size(S,2), size(S,3), M*K);
    fprintf('  Lp = %d,  Ls = %d,  (J=%d, K=%d, M=%d)\n', Lp, Ls, J, K, M);
end


%% ====================================================================
function A = local_load_3d(matpath, which)
%  Load a .mat file and return the single 3-D array inside it, whatever
%  the stored variable is called (auto-detect). Errors clearly if the
%  file holds zero or several candidate arrays.
    s  = load(matpath);
    fn = fieldnames(s);

    % Keep only numeric 3-D arrays as candidates.
    is3d = false(numel(fn),1);
    for i = 1:numel(fn)
        v = s.(fn{i});
        is3d(i) = isnumeric(v) && ndims(v) == 3;
    end
    cand = fn(is3d);

    if isempty(cand)
        % Fall back: a 2-D array (single error-mic case) is acceptable too.
        is2d = false(numel(fn),1);
        for i = 1:numel(fn)
            v = s.(fn{i});
            is2d(i) = isnumeric(v) && ismatrix(v) && ~isscalar(v);
        end
        cand = fn(is2d);
    end

    assert(~isempty(cand), ...
        'load_paths: no numeric array found in %s file "%s".', which, matpath);
    if numel(cand) > 1
        error(['load_paths: %s file "%s" has several arrays (%s). ' ...
               'Set cfg.%s_var to name the one to use.'], ...
               which, matpath, strjoin(cand, ', '), which);
    end
    A = s.(cand{1});
    if ndims(A) == 2          % promote [error x length] to [error x 1 x length]
        A = reshape(A, size(A,1), 1, size(A,2));
    end
end