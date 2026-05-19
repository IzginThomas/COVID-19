function f = rhs_FP(t,Y,x,X_interfaces, para, Nx, dx)
x_interfaces = X_interfaces(2:Nx);
d = para/2 * (1 - x_interfaces.^2).^2; % diffusion term
dprime = -2 * para * x_interfaces .* (1 - x_interfaces.^2); % D', i.e. first derivative of diffusion term
f = zeros(size(Y));
flux = zeros(size(f,1) + 1, size(f,2));
for j = 1:size(Y,2)
    y = Y(:,j);
    B = sum((x_interfaces - x') .* y) * dx;
    lambda = dx *(B + dprime)./d;
    C = (lambda .* d)/dx;
    delta = 1./lambda - 1./(exp(lambda) - 1);
    idx = lambda < 1e-14;
    if ~isempty(idx)
        delta(idx) = 0.5;
    end
    y_tilde = (1 - delta') .* y(2:end) + delta' .* y(1:end-1);
    flux(2:end-1,j) = C' .* y_tilde + d' .* (y(2:end) - y(1:end-1)) ./ dx;
end
f = (flux(2:end,:) - flux(1:end-1,:)) / dx; % Update f with the computed flux values
end