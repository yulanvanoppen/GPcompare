function [p, test] = pval_posterior(m1, m2, S1, S2)
    Dm = m1 - m2;                                       % mean difference
    DS = nearestSPD(S1 + S2);                           % covariance of mean difference
    
    nu = sum(eig(DS) > 0);                              % degrees of freedom
    test = Dm' * (DS \ Dm);                             % effect size
    
    p = chi2cdf(test, nu, 'upper');                     % p-value
end