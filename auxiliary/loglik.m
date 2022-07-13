function p = loglik(mu_par, Sigma_par, t, y, prior)
    warning('off','MATLAB:singularMatrix')          % surpress warnings for singular matrices
    warning('off','MATLAB:eigs:IllConditionedA')

    if isempty(t), p = 0; return, end
    
    N = numel(y) / length(t);                       % number of trajectories

    mu = smoother(mu_par, t);                       % mean
    Sigma = flexible_covariance(t, Sigma_par);      % covariance

    constant = -N/2 * length(t) * log(2*pi);        % denominator terms
    det_term = -N/2 * log(det(Sigma));              % exponentiated term
    exp_term = -1/2 * sum(pagemtimes(reshape(y - mu, [], 1, N), 'transpose', ...
                                     reshape(Sigma \ (y - mu), [], 1, N), 'none'));

                                                    % log likelihood
    p = constant + det_term + exp_term + prior([mu_par Sigma_par]);
    
    warning('on')
end