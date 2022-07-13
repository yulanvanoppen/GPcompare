function [pbeta_mu, pbeta_Sigma] = posterior_beta(t, y, n_meanpar, Sigma_par, prior_sd)

    N = size(y, 2);                                     % number of cells
        
    t = reshape(t, [], 1);                              % basis:
    c = linspace(t(1), t(end), n_meanpar);              % set kernel centers
    h = t(end) / n_meanpar;                             % set kernel width

    X = basis(t, c, h);                                 % parameterized mean and GP covariance + noise
    Sigma = flexible_covariance(t, Sigma_par);

    covar_term = N * X' * (Sigma \ X);                  % compute posterior directly
    mean_term = sum(X' * (Sigma \ y), 2);
                                                        % mean coefficients posterior mean + covariance
    pbeta_Sigma = inv(prior_sd^-2 * eye(n_meanpar) + covar_term);
    pbeta_mu = pbeta_Sigma * mean_term;
end