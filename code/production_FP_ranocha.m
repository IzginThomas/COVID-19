function P = production_FP_ranocha(t,y,x,x_interfaces, para, Nx,dx)
d = para/2 * (1 - x_interfaces.^2).^2; % diffusion term
dprime = -2 * para * x_interfaces .* (1 - x_interfaces.^2); % D', i.e. first derivative of diffusion term
B = sum((x_interfaces - x') .* y) * dx;
lambda = dx *(B + dprime)./d;
C = (lambda .* d)/dx;
delta = 1./lambda - 1./(exp(lambda) - 1);
idx = lambda < 1e-14;
if ~isempty(idx)
    delta(idx) = 0.5;
end
y_tilde = (1 - delta) .* y + delta .* circshift(y,1);
P = sparse(Nx,Nx);
for i = 2:(Nx-1)
    P(i, i + 1) = (max(0, C(i + 1)) * y_tilde(i+1,i + 1) + d(i + 1) * y(i + 1)/dx)/dx;
    P(i, i - 1) = (-min(0, C(i)) * y_tilde(i,i) + d(i) * y(i - 1)/dx)/dx;
end
P(1, 2) = (max(0, C(2)) * y_tilde(2,2) + d(2) * y(2)/dx)/dx;
P(Nx, Nx - 1) = (-min(0, C(Nx)) * y_tilde(Nx,Nx) + d(Nx) * y(Nx - 1)/dx)/dx;

end