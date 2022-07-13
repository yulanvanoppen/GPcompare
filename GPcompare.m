% GPcompare Compare groups of time series using Guassian processes (GPs) by
% evaluating their posterior mean difference distribution.
%
% [p, T, E] = GPcompare(X, Y, ...) compares two groups of time series data
% given in 2-D arrays X and Y with rows corresponding to (aligned) time
% points and columns to individuals (or repetitions). The dimensions of X 
% and Y are not restricted to be equal.
%
% Output:
% - [p] tail probability of the zero vector;
% - [T] effect size;
% - [E] effect size decomposition such that T = E'E.
%
function [p, T, E] = GPcompare(X, Y, varargin)
    addpath('auxiliary')

    %% Parse input
    defaultYLimits = [min(min(X, Y), [], 'all') max(max(X, Y), [], 'all')];
    defaultLim = [];
    defaultPlot = true;
    
    checkData = @(data) isnumeric(data) && ~isempty(data);
    checkLimit = @(lim) isnumeric(lim) && (isempty(lim) || numel(lim) == 1 && lim > 0);
    checkYLimits = @(lim) isnumeric(lim) && numel(lim) == 2 && lim(1) < lim(2);
    
    parser = inputParser;
    addRequired(parser, 'X', checkData)
    addRequired(parser, 'Y', checkData)
    addOptional(parser, 'plot', defaultPlot, @islogical)
    addOptional(parser, 'yLimits', defaultYLimits, checkYLimits)
    addOptional(parser, 'differenceLimit', defaultLim, checkLimit)
    addOptional(parser, 'effectLimit', defaultLim, checkLimit)
    
    parser.KeepUnmatched = true;
    parse(parser, X, Y, varargin{:})
    

    %% Settings
    n_meanpar = 10;                                     % numbers of basis functions
    n_Sigmapar = 5;

    tX = linspace(0, 1, size(X, 1));                    % data time grids
    tY = linspace(0, 1, size(Y, 1));
    t_coarse = linspace(0, 1, n_meanpar);               % time grid for comparisons
    t_fine = linspace(0, 1, 100);                       % time grid for plotting


    %% Fitting
                                                        % coefficients ~ N(0, 10^2)
    prior = @(par) sum(log(normpdf(par(1:n_meanpar), 0, 10))) ...     
                    ...                                 % length_scales L(t) coefficients ~ Gamma(1/5, 1)
                    + sum(log(gampdf(exp(par(n_meanpar+n_Sigmapar+1:end-1)), 1/5, 1))) ...
                    ...                                 % amplitudes and noise coefficients ~ G(1, 1)
                    + sum(log(gampdf(exp(par([n_meanpar+1:n_meanpar+n_Sigmapar, end])), 1, 1)));

    [optX, ~, ~] = fit(tX, X, prior, n_meanpar, n_Sigmapar);
    [optY, ~, ~] = fit(tY, Y, prior, n_meanpar, n_Sigmapar);

                                                        % compute mean coefficients' posterior
    [pbeta_mX, pbeta_SX] = posterior_beta(tX, X, n_meanpar, optX(n_meanpar+1:end), 10);
    [pbeta_mY, pbeta_SY] = posterior_beta(tY, Y, n_meanpar, optY(n_meanpar+1:end), 10);

                                                        % transform posterior beta to GP mean
    [pmu_mX, pmu_SX] = posterior_mu(t_coarse, pbeta_mX, pbeta_SX);
    [pmu_mY, pmu_SY] = posterior_mu(t_coarse, pbeta_mY, pbeta_SY);
    [pmu_mX_fine, pmu_SX_fine] = posterior_mu(t_fine, pbeta_mX, pbeta_SX);
    [pmu_mY_fine, pmu_SY_fine] = posterior_mu(t_fine, pbeta_mY, pbeta_SY);

    
    %% Compare
    pDmu_m = pmu_mX - pmu_mY;                           % posterior mean difference (mean + covariance)
    pDmu_S = nearestSPD(pmu_SX + pmu_SY);
    E = chol(pDmu_S)' \ pDmu_m;                         % effect size decomposition
    
                                                        % tail probability and effect size
    [p, T] = pval_posterior(pmu_mX, pmu_mY, pmu_SX, pmu_SY);

    pDmu_m_fine = pmu_mX_fine - pmu_mY_fine;            % posterior mean difference (fine grid)
    pDmu_S_fine = nearestSPD(pmu_SX_fine + pmu_SY_fine);
    pDmu_sd_fine = sqrt(diag(pDmu_S_fine));
    
    
    %% Plotting
    if parser.Results.plot
        if isempty(parser.Results.differenceLimit)          % set y-axis limits for panel 3
            differenceLimit = max(0.000001, max(max(abs(pDmu_m_fine + 2*pDmu_sd_fine)), ...
                                                max(abs(pDmu_m_fine - 2*pDmu_sd_fine))));
        else
            differenceLimit = parser.Results.differenceLimit;
        end
        if isempty(parser.Results.effectLimit)
            effectLimit = max(0.000001, max(abs(E)));
        else
            effectLimit = parser.Results.effectLimit;
        end


        f = figure;                                         % initialize 3-panel figure
        f.Position = [100 100 1000 400];
        tiledlayout(1, 3)


        nexttile(1)                                         % first data set + mean posterior
        plot(tX, X, 'Color', [0.9290, 0.6940, 0.1250, .5]);
        hold on
        plot(t_fine, pmu_mX_fine, 'LineWidth', 2, 'Color', [0.9290, 0.6940, 0.1250])
        plot(t_fine, [pmu_mX_fine + 2*sqrt(diag(pmu_SX_fine)), pmu_mX_fine - 2*sqrt(diag(pmu_SX_fine))], ...
             '--', 'LineWidth', 2, 'Color', [0.9290, 0.6940, 0.1250])
        patch([t_fine fliplr(t_fine)], ...
              [(pmu_mX_fine + 2*sqrt(diag(pmu_SX_fine)))' fliplr((pmu_mX_fine - 2*sqrt(diag(pmu_SX_fine)))')], ...
              'g', 'FaceColor', [0.9290, 0.6940, 0.1250], 'FaceAlpha', .3, 'EdgeAlpha', 0)      
        ylim(parser.Results.yLimits)
        ylabel('measurement')
        xlabel('cell cycle progression')
        xticks(0:.2:1)
        title('X')


        nexttile(2)                                         % second data set + mean posterior
        plot(tY, Y, 'Color', [0.4940, 0.1840, 0.5560, .5]);
        hold on
        plot(t_fine, pmu_mY_fine, 'LineWidth', 2, 'Color', [0.4940, 0.1840, 0.5560])
        plot(t_fine, [pmu_mY_fine + 2*sqrt(diag(pmu_SY_fine)), pmu_mY_fine - 2*sqrt(diag(pmu_SY_fine))], ...
             '--', 'LineWidth', 2, 'Color', [0.4940, 0.1840, 0.5560])
        patch([t_fine fliplr(t_fine)], ...
              [(pmu_mY_fine + 2*sqrt(diag(pmu_SY_fine)))' fliplr((pmu_mY_fine - 2*sqrt(diag(pmu_SY_fine)))')], ...
              'g', 'FaceColor', [0.4940, 0.1840, 0.5560], 'FaceAlpha', .3, 'EdgeAlpha', 0)    
        ylim(parser.Results.yLimits)
        ylabel('measurement')
        xlabel('cell cycle progression')
        xticks(0:.2:1)
        title('Y')


        nexttile(3)                                         % posterior mean difference and effect size decomposition
        yyaxis left
        plot(t_fine, pDmu_m_fine, 'LineWidth', 2, 'Color', [0 0.4470 0.7410])
        hold on
        plot(t_fine, [pDmu_m_fine + 2*pDmu_sd_fine, pDmu_m_fine - 2*pDmu_sd_fine], ...
             '--', 'LineWidth', 2, 'Color', [0 0.4470 0.7410])
        patch([t_fine fliplr(t_fine)], [(pDmu_m_fine + 2*pDmu_sd_fine)' fliplr((pDmu_m_fine - 2*pDmu_sd_fine)')], ...
              'g', 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', .3, 'EdgeAlpha', 0)
        ylim([-1 1] * differenceLimit * 1.05)
        ylabel('measurement difference')
        xlabel('cell cycle progression')
        xticks(0:.2:1)

        yyaxis right
        bar(t_coarse, E, 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', .3, 'EdgeColor', [0.8500 0.3250 0.0980])
        ylim([-1 1] * effectLimit * 1.05)
        ylabel('effect size contribution')

        title('X - Y')
    end
end