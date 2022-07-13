function [par, opt, H] = fit(t, y, prior, n_meanpar, n_Sigmapar)
    opt = Inf;
    n_covpar = n_Sigmapar*2 + 1;
    par = zeros(1, n_meanpar+n_covpar);                 % initialize maximizer and its Hessian matrix
    H = zeros(n_meanpar+n_covpar);
    
    options = optimoptions('fmincon', 'Display', 'off', 'OptimalityTolerance', 1e-4);
    for n = 1:5
        try                                             % maximize likelihood a posteriori
            [MLE_par, MLE, ~, ~, ...
             ~, ~, hess] = fmincon(@(par) -loglik(par(1:n_meanpar), par(n_meanpar+1:end), t, y, prior), ...
                                   randn(1, length(par))+1, [], [], [], [], ...     % initial value
                                   -5*ones(1, length(par)), ...                     % lbound
                                   5*ones(1, length(par)), [], options);            % ubound
        catch
            continue;
        end
        if MLE<opt                                      % save on improvement
            par = MLE_par;
            opt = MLE;
            H = hess;
        end
    end
end