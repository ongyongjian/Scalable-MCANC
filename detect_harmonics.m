function [f0, harm_freqs] = detect_harmonics(cfg, x, varargin)
%% detect_harmonics  Auto-detect fundamental frequency and harmonics
%
%  WHAT IT DOES
%  ------------
%  Finds the strongest spectral peak within [cfg.f_lo, cfg.f_hi] and treats
%  it as the fundamental f0. Returns f0 and its integer multiples up to
%  cfg.f_hi, for marking harmonics on PSD/spectrum plots.
%
%  General J x K x M: x may be [J x N] (multiple reference sensors). The
%  fundamental is a property of the common source, so by default the peak
%  is found on the AVERAGED reference spectrum across all J sensors, which
%  is more robust than any single sensor. Set 'ref' to use one specific
%  sensor instead. Reduces to the old behaviour when J = 1.
%
%  INPUTS
%  ------
%  cfg      : struct — needs fs, f_lo, f_hi
%  x        : [J x N]  reference signal(s)
%  varargin : optional name-value pairs
%             'f0'      , scalar  — supply f0 directly to skip detection
%             'n_harm'  , integer — max number of harmonics (default 8)
%             'ref'     , integer — use only this reference row for detection
%                                  (default 0 = average all J references)
%
%  OUTPUTS
%  -------
%  f0         : scalar    detected (or supplied) fundamental frequency (Hz)
%  harm_freqs : [1 x Nh]  harmonic frequencies f0, 2*f0, ... up to cfg.f_hi
%
%  CALL EXAMPLE
%  ------------
%  [f0, harm_freqs] = detect_harmonics(cfg, x);
%  [f0, harm_freqs] = detect_harmonics(cfg, x, 'ref', 1);

p = inputParser;
addParameter(p, 'f0',     [], @isnumeric);
addParameter(p, 'n_harm', 8,  @isnumeric);
addParameter(p, 'ref',    0,  @isnumeric);
parse(p, varargin{:});
f0     = p.Results.f0;
n_harm = p.Results.n_harm;
ref    = p.Results.ref;

J = size(x, 1);

if isempty(f0)
    if ref >= 1 && ref <= J
        % Use a single chosen reference sensor.
        [Pxx, F] = pwelch(x(ref,:), hann(8192), 4096, 8192, cfg.fs);
    else
        % Average the PSD across all J reference sensors (default).
        Pxx = 0;  F = [];
        for j = 1:J
            [Pj, F] = pwelch(x(j,:), hann(8192), 4096, 8192, cfg.fs);
            Pxx = Pxx + Pj;
        end
        Pxx = Pxx / J;
    end
    band_mask = F >= cfg.f_lo & F <= cfg.f_hi;
    [~, idx0] = max(Pxx(band_mask));
    F_band    = F(band_mask);
    f0        = F_band(idx0);
end

harm_freqs = f0 * (1:n_harm);
harm_freqs = harm_freqs(harm_freqs <= cfg.f_hi);
end