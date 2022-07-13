function C = flexible_covariance(t, par)
    n_coef = (length(par)-1) / 2;
    if n_coef > 1
        amplitudes = smoother(exp(par(1:n_coef)), t);   % variable GP amplitude
        coefs2 = exp(par(n_coef+1:end-1));
        length_scales = smoother(coefs2, t);            % variable length scale
    else
        amplitudes = exp(par(1));
        length_scales = exp(par(2));
    end
    noise = exp(par(end));
    C = covariance_function(t - t', amplitudes, length_scales, noise);
end

function K = covariance_function(t, amplitudes, length_scales, noise)
    nu = 5/2;   
    denom = sqrt((length_scales.^2 + length_scales'.^2)/2);
    inner = sqrt(2*nu) * abs(t - t') ./ denom;
                                                        % variable length scale normalization
    normalizer = sqrt(length_scales .* length_scales') ./ denom;
    
                                                        % matern correlation function
    K = normalizer .* 2^(1-nu) / gamma(nu) .* inner.^nu .* besselk(nu, inner);
    K(logical(eye(size(K)))) = 1;
    
    K(K < 1e-12) = 0;                                   % nullify negligible elements
    K = amplitudes * amplitudes' .* K + noise^2 * eye(size(K));   % scale according to variance
end