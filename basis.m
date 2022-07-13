function b = basis(x, c, h)
    b = exp(-1/2 * (x-c).^2 / h^2);                     % normal kernels
end