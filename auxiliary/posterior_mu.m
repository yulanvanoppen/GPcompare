function [pmu_mean, pmu_Sigma] = posterior_mu(t, pbeta_mu, pbeta_Sigma)
    t = reshape(t, [], 1);                              % basis:
    c = linspace(t(1), t(end), length(pbeta_mu));       % set kernel centers
    h = t(end) / length(pbeta_mu);                      % set kernel width

    X = basis(t, c, h);                                 % basis vectors
    
    pmu_mean = X * pbeta_mu;                            % linear combination of beta's posterior
    pmu_Sigma = X * pbeta_Sigma * X';
end