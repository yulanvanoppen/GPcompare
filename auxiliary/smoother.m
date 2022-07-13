function y = smoother(coefs, t)
    coefs = reshape(coefs, 1, []);                      % basis coefficients
    
    t = reshape(t, [], 1);                              % force to column vector
    c = linspace(t(1), t(end), length(coefs));          % set kernel centers
    h = t(end) / length(coefs);                         % set kernel width
    
    y = sum(coefs .* basis(t, c, h), 2);                % linear combination of basis functions
end