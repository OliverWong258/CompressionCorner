function lamda_y = MaxEigY(U,gamma)

rho = U(:,:,1);
u = U(:,:,2)./rho;
v = U(:,:,3)./rho;
E = U(:,:,4);
p = (gamma-1)*(E - 0.5*rho.*(u.^2+v.^2));
c = sqrt(gamma*p./rho);
lamda_y = max( max(abs([v;v-c; v+c]),[],2) );

end

