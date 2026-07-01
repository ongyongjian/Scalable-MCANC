function [d_s, e_s, n_valid] = nan_safe(d, e)
%% nan_safe  Trim d and e to the columns where BOTH are finite, so spectral
%            functions (pwelch, spectrogram, cpsd) never see NaN/Inf after a
%            divergence. Returns the trimmed signals and the valid length.
%
%  If the controller diverged at sample n_div, e(:,n_div:end) is NaN. This
%  keeps only the leading finite block so the plots still work (showing the
%  pre-divergence behaviour) instead of crashing.
%
%  d, e : [M x N]   ->   d_s, e_s : [M x n_valid]

    finite_cols = all(isfinite(d), 1) & all(isfinite(e), 1);
    if all(finite_cols)
        d_s = d;  e_s = e;  n_valid = size(d,2);
        return;
    end
    % first non-finite column = end of the usable block
    first_bad = find(~finite_cols, 1, 'first');
    if isempty(first_bad), first_bad = size(d,2)+1; end
    n_valid = first_bad - 1;
    if n_valid < 1
        error('nan_safe: signals are non-finite from the very first sample.');
    end
    d_s = d(:, 1:n_valid);
    e_s = e(:, 1:n_valid);
    warning('nan_safe:trimmed', ...
        ['Divergence detected: using first %d of %d samples for this plot ' ...
         '(pre-divergence window).'], n_valid, size(d,2));
end