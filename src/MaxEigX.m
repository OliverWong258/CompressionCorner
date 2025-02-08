function lamda_x = MaxEigX(U,gamma)
rho = U(:,:,1);
u = U(:,:,2)./rho;
v = U(:,:,3)./rho;
E = U(:,:,4);
p = (gamma-1)*(E - 0.5*rho.*(u.^2+v.^2));
c = sqrt(gamma*p./rho);
lamda_x = max( max(abs([u;u-c; u+c]),[],2) );

end