function U = genLaguerreGaussian(X,Y,z,wavlen,w0,l,p)

n = 1;
zR = pi*w0^2*n/wavlen;
N = abs(l) + 2*p;
k = 2*pi*n/wavlen;

w = @(z) w0*sqrt(1 + (z/zR).^2);
psi = @(z) (N+1)*atan(z/zR);

[phi,rho] = cart2pol(X,Y);

% C = sqrt(2*factorial(p)/(pi*factorial(p+abs(l))));

C = 1;
U = C .* 1./w(z) * (rho*sqrt(2)./w(z)).^abs(l) .*exp(-rho.^2./w(z).^2) ...
    .* myLaguerre(p,abs(l),2*rho.^2./w(z).^2)...
    .*exp(-1i*k*rho.^2./2./R(z)).*exp(-1i*l*phi).*exp(1i*psi(z));

U = U/sqrt(sum(abs(U(:)).^2));

function v = R(z)
if z == 0
    v = inf;
else
    v = z*(1 + (zR./z).^2);
end
end

end